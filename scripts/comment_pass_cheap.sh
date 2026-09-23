#!/usr/bin/env bash
# =============================================================================
# Passe de commentaires ÉCONOMIQUE et automatique.
#
# Principe : un seul appel à Claude par fichier (modèle Sonnet, aucun outil, un seul
# tour). Le script lui envoie la charte + la liste des red flags + le fichier, et Claude
# renvoie le fichier réécrit. Le script vérifie lui-même que seul du commentaire a
# changé, puis écrit le fichier et commite. Pas d'agent qui relit le projet en boucle.
#
# Usage (à la racine du dépôt, sur ta branche de travail, arbre propre) :
#     bash scripts/comment_pass_cheap.sh                       # dossiers par défaut
#     bash scripts/comment_pass_cheap.sh src/solver/fe/impl/common   # dossiers choisis
#
# Relançable : les fichiers déjà commités sont sautés.
# Variables : MODEL (défaut sonnet), MAX_LINES (défaut 1500, fichiers plus longs sautés).
# =============================================================================
set -uo pipefail

die() { echo; echo "!! $*"; exit 1; }

cd "$(git rev-parse --show-toplevel 2>/dev/null)" || die "Pas dans un dépôt git."
ROOT=$(pwd)
MODEL=${MODEL:-sonnet}
MAX_LINES=${MAX_LINES:-1500}
LOG=logs/comment-pass-cheap
mkdir -p "$LOG" docs
touch docs/design-red-flags.md
grep -qxF "logs/" .git/info/exclude 2>/dev/null || echo "logs/" >> .git/info/exclude

DIRS=("$@")
if [ ${#DIRS[@]} -eq 0 ]; then
  DIRS=(src/model/mesh/api src/solver/fe/api src/gradient/api
        src/discretization/fe/api src/io/api src/utils)
fi

# ------------------------------------------------------------------ vérifications
for t in claude git gcc awk; do command -v "$t" >/dev/null || die "Outil manquant : $t"; done
case "$(git rev-parse --abbrev-ref HEAD)" in
  main|master) die "Tu es sur main : passe sur ta branche de travail." ;;
esac
if ! git diff --quiet || ! git diff --cached --quiet; then
  die "Des fichiers suivis sont modifiés. Commite-les ou annule-les (git status) puis relance."
fi

CHARTER=$(awk '/^## Comments and documentation/{p=1; print; next} p && /^## /{exit} p' CONTRIBUTING.md)
[ -n "$CHARTER" ] || die "Section '## Comments and documentation' introuvable dans CONTRIBUTING.md."

# Code sans commentaires ni espaces : deux fichiers qui donnent la même chose ne
# diffèrent que par leurs commentaires.
strip() { gcc -fpreprocessed -dD -E -P -x c++ - 2>/dev/null | tr -d '[:space:]'; }

# Pas de grep -q : avec pipefail, sa sortie anticipée ferait échouer git log (SIGPIPE).
already_done() { git log --format=%s | grep -xF "$1" > /dev/null; }

# ------------------------------------------------------------------ étape 0 : ASCII, sans IA
if ! already_done "comments: ASCII only"; then
  echo "== Étape 0 : remplacement des caractères Unicode (sans IA)"
  git ls-files -z 'src/*.h' 'src/*.hpp' 'src/*.cc' 'src/*.cpp' 'src/*.cpp.in' | \
    xargs -0 -r sed -i 's/→/->/g; s/←/<-/g; s/⇒/=>/g; s/—/ - /g; s/–/-/g; s/≤/<=/g; s/≥/>=/g; s/≠/!=/g; s/×/x/g; s/²/^2/g; s/³/^3/g; s/·/*/g'
  bad=0
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    cmp -s <(git show "HEAD:$f" | strip) <(strip < "$f") || { echo "   caractère dans du code : $f"; bad=1; }
  done < <(git diff --name-only)
  if [ $bad -eq 1 ]; then
    git checkout -- src
    echo "   Remplacement annulé (un caractère était dans du code, pas dans un commentaire)."
  else
    git add -A src
    git commit -q --allow-empty -m "comments: ASCII only"
    echo "   OK, commité."
  fi
fi

# ------------------------------------------------------------------ consignes envoyées à Claude
read -r -d '' INSTRUCTIONS <<'EOF'
You edit the comments of ONE C++ source file of FUnTiDES, a C++17/Kokkos
spectral-element wave-propagation solver. Do not use any tool: everything you need
is in this message.

Task: apply the charter below to this file. Keep, rewrite or delete each comment.
Add the missing Doxygen interface comments on public entities. Make the Doxygen
format uniform. Do not add references to docs/design.md.

ABSOLUTE RULE: change nothing outside comments. Not one character of code, no
reformatting, no reordering, no include change. Your output is compared with the
original after stripping all comments and whitespace; it is rejected if anything
differs.

Never guess: if a meaning, unit or convention is uncertain, write
`@todo VERIFY: <precise question>`.
Suspected bugs and design problems: do not fix them and do not hide them in a
comment; list them in the REDFLAGS block. Do not repeat the known red flags.

Reply with exactly this format and nothing else:
<<<FILE
(the complete new content of the file)
FILE>>>
<<<REDFLAGS
- `path`, `symbol`: problem in one or two sentences.
(leave empty if there is none)
REDFLAGS>>>
EOF

# ------------------------------------------------------------------ passe fichier par fichier
FILES=()
for d in "${DIRS[@]}"; do
  [ -d "$d" ] || { echo "-- $d : dossier absent, sauté"; continue; }
  while IFS= read -r f; do FILES+=("$f"); done < <(
    git ls-files "$d" | grep -E '\.(h|hpp|cc|cpp|cpp\.in)$' | grep -v 'LagrangeBasis')
done

todo=0; total_lines=0
for f in "${FILES[@]}"; do
  already_done "comments: $f" && continue
  todo=$((todo + 1)); total_lines=$((total_lines + $(wc -l < "$f")))
done
echo
echo "== $todo fichiers à traiter, $total_lines lignes au total, modèle $MODEL"

done_count=0
for f in "${FILES[@]}"; do
  already_done "comments: $f" && continue
  n=$(wc -l < "$f")
  if [ "$n" -gt "$MAX_LINES" ]; then
    echo "-- $f : $n lignes (> $MAX_LINES), sauté" | tee -a "$LOG/skipped.txt"
    continue
  fi

  done_count=$((done_count + 1))
  echo "== [$done_count/$todo] $f ($n lignes)"
  tag=${f//\//_}
  prompt="$LOG/$tag.prompt"; out="$LOG/$tag.out"; new="$LOG/$tag.new"

  {
    printf '%s\n\n' "$INSTRUCTIONS"
    printf '## Charter\n%s\n\n' "$CHARTER"
    printf '## Known red flags\n'; cat docs/design-red-flags.md
    printf '\n## File: %s\n' "$f"
    cat "$f"
    printf '\n'
  } > "$prompt"

  # Lancé depuis /tmp pour ne pas charger CLAUDE.md ni les réglages du projet :
  # Claude ne reçoit que ce que le script lui envoie.
  if ! (cd /tmp && claude -p "Follow the instructions given in the input." \
          --model "$MODEL" --max-turns 1 --output-format text \
          --disallowedTools Bash Read Edit Write Glob Grep WebFetch WebSearch Task Agent \
          < "$ROOT/$prompt") > "$out" 2> "$LOG/$tag.err"; then
    if grep -qi "max turns" "$out" "$LOG/$tag.err"; then
      echo "   a tenté d'utiliser un outil, sauté"; echo "$f : a tenté d'utiliser un outil" >> "$LOG/skipped.txt"; continue
    fi
    cat "$LOG/$tag.err" "$out" | tail -5
    die "Claude s'est arrêté sur $f (limite d'usage probable). Relance le script plus tard : il reprendra ici."
  fi

  if ! grep -qx '<<<FILE' "$out" || ! grep -qx 'FILE>>>' "$out"; then
    echo "   réponse incomplète (fichier trop long ?), sauté"; echo "$f : réponse incomplète" >> "$LOG/skipped.txt"
    continue
  fi
  awk '/^<<<FILE$/{p=1; next} /^FILE>>>$/{p=0} p' "$out" > "$new"

  if ! cmp -s <(strip < "$f") <(strip < "$new"); then
    echo "   Claude a modifié du code : rejeté, fichier inchangé"; echo "$f : code modifié, rejeté" >> "$LOG/skipped.txt"
    continue
  fi

  cp "$new" "$f"
  awk '/^<<<REDFLAGS$/{p=1; next} /^REDFLAGS>>>$/{p=0} p' "$out" | grep -E '^\s*-' >> docs/design-red-flags.md
  git add "$f" docs/design-red-flags.md
  git commit -q --allow-empty -m "comments: $f"
  echo "   OK, commité ($(grep -c '@todo VERIFY' "$f") @todo VERIFY dans le fichier)"
done

# ------------------------------------------------------------------ bilan
echo
echo "== Terminé"
echo "Fichiers sautés (à voir à la main) : $LOG/skipped.txt"
echo "Points à vérifier                  : grep -rn '@todo VERIFY' src"
echo "Red flags                          : docs/design-red-flags.md"
echo "Commits                            : git log --oneline main..HEAD"
