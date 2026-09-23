# CLAUDE.md

@CONTRIBUTING.md

## Current task: comment and documentation pass

- Only comments and documentation change. Code, formatting, names and includes stay
  byte-for-byte identical outside comments.
- `scripts/check_comments_only.sh` proves it: it compares the code with comments
  stripped by gcc, so it also catches any comment edit that would change what the
  compiler sees. It must pass before any commit. Do not build the project: the check
  is sufficient and a build is very slow.
- Doxygen (`doxygen Doxyfile`, warnings in `build/doxygen-warnings.log`) runs once, at
  the end of the whole pass.
- Cross-module conventions go in `docs/design.md`; design problems and suspected bugs
  go in `docs/design-red-flags.md` (read it first: some are already known).
