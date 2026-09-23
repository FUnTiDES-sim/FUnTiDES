Comment pass on the directories given here: $ARGUMENTS
If none are given, use: src/model/mesh/api src/solver/fe/api src/gradient/api
src/discretization/fe/api src/io/api src/utils

You orchestrate only: never read or edit source files yourself, do not build, do not
ask me questions, keep your messages to one line per file.

1. Run once: `git log --format=%s main..HEAD` (a file F is done if a commit is titled
   "comments: F").
2. Run once:
   `git ls-files <directories> | grep -E '\.(h|hpp|cc|cpp|cpp\.in)$' | grep -v LagrangeBasis | xargs wc -l`
   Keep the files that are not done and have at most 1500 lines. The others are skipped.
3. For each kept file F, in order:
   a. Delegate F to the comment-file subagent, passing only the path.
   b. Run this single command (replace F):
      `scripts/check_comments_only.sh && git add F docs/design-red-flags.md && git commit -q -m "comments: F" && echo DONE || { git checkout -- F; echo REJECTED; }`
   c. Print one line: F, then DONE or REJECTED.
4. At the end, print the skipped files and the rejected files.
