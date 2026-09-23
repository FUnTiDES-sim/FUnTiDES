#!/usr/bin/env bash
# Fails if a C++ file changed since <ref> differs in anything other than comments.
# Default <ref>: the merge-base with main.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
ref=${1:-$(git merge-base main HEAD)}
strip() { gcc -fpreprocessed -dD -E -P -x c++ - 2>/dev/null | tr -d '[:space:]'; }
status=0
while IFS= read -r f; do
  [ -n "$f" ] || continue
  if ! cmp -s <(git show "$ref:$f" | strip) <(strip < "$f"); then
    echo "CODE CHANGED: $f"; status=1
  fi
done < <(git diff --name-only --diff-filter=M "$ref" -- '*.h' '*.hpp' '*.cc' '*.cpp' '*.cpp.in')
while IFS= read -r f; do
  [ -n "$f" ] || continue
  echo "C++ FILE DELETED OR RENAMED: $f"; status=1
done < <(git diff --name-only --diff-filter=DR "$ref" -- '*.h' '*.hpp' '*.cc' '*.cpp' '*.cpp.in')
[ $status -eq 0 ] && echo "OK: only comments changed since ${ref:0:10}"
exit $status
