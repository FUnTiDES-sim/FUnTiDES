# Contributing Guidelines

This document defines the mandatory code formatting standards and submission requirements for this project. All contributions must conform to the specified formatting rules to be accepted.

## Code Formatting

This project uses **clang-format** to enforce consistent code style. All C++ code must be formatted according to our style configuration before submission.

### Formatting Style

Our code style is based on **Google C++ Style Guide** with the following custom modifications:

- **Brace Style**: Custom brace wrapping
  - Braces are placed on new lines for classes, functions, structs, unions, namespaces, enums
  - Braces are placed on new lines for control statements (if, for, while, etc.)
  - Braces are placed on new lines before catch and else statements
- **Indentation**: 2 spaces (no tabs)
- **Column Limit**: 80 characters per line

### Before Submitting a Pull Request

**1. Install clang-format**

We use version 18 of clang-format.

```bash
# Ubuntu/Debian
sudo apt install clang-format-18

# macOS
brew install clang-format@18

# Python (cross-platform)
#   - clang only 
pip install clang-format==18
#   - all dev requirements at once
pip install -r requirements-dev.txt

# Or check if already installed
clang-format --version
```

**2. Format your code**

Before committing any changes, run clang-format on all modified C++ files:

```bash
# Format a single file
clang-format -i path/to/your/file.cpp

# Format all C++ files recursively
find . \( -name "*.h" -o -name "*.hpp" -o -name "*.cc" -o -name "*.cpp" -o -name "*.cxx" \) -exec clang-format -i {} \;
```

**3. Verify formatting**

Check that your files are properly formatted:

```bash
# Check if files need formatting (should show no output if properly formatted)
find . \( -name "*.h" -o -name "*.hpp" -o -name "*.cc" -o -name "*.cpp" -o -name "*.cxx" \) -exec clang-format -i {} \;
git diff --exit-code
```

### Configuration Details

The project's code formatting is defined in [`.clang-format`](.clang-format):

## CI/CD

### Format Checking

Our continuous integration system automatically checks code formatting.
**Pull requests with formatting violations will fail the build and stop any subsequent build**. 

If you see formatting errors in CI:

1. Pull the latest changes from your branch
2. Run clang-format as described above
3. Commit and push the formatting fixes

### Caching

Our continuous integration system uses sccache to speed up the builds.
In some rare instances where the cache might be corrupted, you can deactivate it by adding the tag `SCCACHE_OFF` to your PR.

## Submission Checklist

Before submitting your pull request:

- [ ] Code is properly formatted with clang-format
- [ ] All tests pass locally
- [ ] Commit messages are clear and descriptive
- [ ] No formatting violations in CI checks

## Questions?

If you encounter any issues with formatting or have questions about the style guide, please open an issue or ask in the pull request discussion.

Thank you for helping maintain code quality!
## Comments and documentation

Principles (after Ousterhout, *A Philosophy of Software Design*):

1. A comment says what the code cannot say: why, invariants, units, index
   conventions, preconditions, ownership, host/device constraints.
   A comment that paraphrases the next line is deleted.
2. Interface comments (headers, classes, public functions and members) describe
   what the entity does and how to use it, never how it is implemented.
   Every public entity under */api/ has one.
3. Implementation comments (function bodies) are rare: they explain the intent of
   a non-obvious block, not each line.
4. Low-level comments add precision (units, sizes, layouts: "size numNodes,
   indexed by linearIndex(i,j,k)"); high-level comments add intuition (role
   of the class in the time loop).
5. Cross-module design decisions (index conventions, GPU dispatch rules,
   explicit instantiation scheme) live once in docs/design.md; code refers
   to them with @see instead of repeating them.
6. No history in comments ("ported from", "V1", "for retrocompatibility"):
   that belongs in git log. A TODO states what is missing and why.
7. If a clear interface comment cannot be written in about three sentences, the
   design is probably wrong: do not paper over it, log it in
   docs/design-red-flags.md.
8. Never guess semantics. If the meaning of a parameter or a unit is uncertain,
   write `@todo VERIFY: <question>` and move on.
9. Interface comments never name the implementations (Doxygen lists derived
   classes) nor describe what a particular implementation or caller does. If the
   behavior genuinely differs between implementations, log it as a red flag.
10. A suspected bug (wrong type, swapped template arguments, dead parameter) is a
    red flag entry, not a @todo VERIFY, and its comment must not hide it.
11. A comment that describes an architecture not yet implemented (unused base
    class, unwired traits) is replaced by a factual note and a red flag entry.

Format:
- Doxygen `/** ... */` blocks with @brief for classes and functions,
  `///<` for members; @tparam, @param[in|out|in,out], @return, @throws.
- A @param that only restates the parameter name is omitted.
- English, ASCII only (no arrows, dashes or math symbols), no banner lines.
