Run the comment and documentation pass on the whole code base. You orchestrate; you do
not write comments yourself. Do not build the project at any point: the check script is
sufficient. Do not ask me questions.

Directories, in this order:
src/model/mesh/api, src/solver/fe/api, src/gradient/api, src/discretization/fe/api,
src/io/api, src/utils, src/discretization/fe/impl/common,
src/discretization/fe/impl/makutu, src/discretization/fe/impl/tensorial,
src/model/mesh/impl, src/solver/fe/impl/common, src/solver/fe/impl/acoustic,
src/solver/fe/impl/elastic, src/solver/fe/impl/acoustoelastic, src/solver/fe/DG,
src/solver/fe/DG-SEm, src/solver/fe/DG_padaptive, src/gradient/impl, src/io/impl,
src/main, src/model/mesh/pywrap, src/solver/fe/pywrap, src/gradient/pywrap

For each directory D that exists and has no commit titled "comments: D" in
`git log main..HEAD`:
1. Delegate D to the comment-writer subagent (one subagent call per directory).
2. Run `scripts/check_comments_only.sh`. If it fails, send the output back to a new
   comment-writer call on D to fix it, once. If it still fails, run
   `git checkout -- D docs`, record D as skipped, and go to the next directory.
3. `git add D docs` then `git commit -m "comments: D"`.
4. Print one line: D, number of @todo VERIFY added, number of red flags added.

At the end: run `mkdir -p build && doxygen Doxyfile`, then print the number of lines of
build/doxygen-warnings.log, the total `grep -rn "@todo VERIFY" src | wc -l`, the
skipped directories, and `git log --oneline main..HEAD`.
