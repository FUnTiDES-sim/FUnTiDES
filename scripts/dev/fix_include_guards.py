#!/usr/bin/env python3
"""
Script to fix include guards according to Google C++ style guide.

Google standard for include guards:
- Format: <PROJECT>_<PATH>_<FILE>_H_
- All uppercase
- Based on the full path in the project (from src/)
- Ends with _H_

Example: foo/src/bar/baz.h in project foo -> FOO_BAR_BAZ_H_
"""

import os
import re
import sys
from pathlib import Path
from typing import Tuple, Optional


def should_exclude_path(filepath: Path, src_root: Path, exclude_dirs: list) -> bool:
    """
    Check if a file should be excluded based on its path.
    
    Args:
        filepath: Path of the file to check
        src_root: Root of the src folder
        exclude_dirs: List of directory names to exclude
    
    Returns:
        True if the file should be excluded, False otherwise
    """
    if not exclude_dirs:
        return False
    
    try:
        rel_path = filepath.relative_to(src_root)
    except ValueError:
        return False
    
    # Check if any parent directory is in the exclusion list
    for part in rel_path.parts[:-1]:  # Exclude the filename
        if part in exclude_dirs:
            return True
    
    return False


def generate_guard_name(filepath: Path, src_root: Path, project_name: str = "") -> str:
    """
    Generate the guard name according to Google standard.
    
    Args:
        filepath: Path of the header file
        src_root: Root of the src folder
        project_name: Project name (optional)
    
    Returns:
        The guard name (e.g., FUNTIDES_SOLVER_FE_MESH_H_)
    """
    # Get relative path from src/
    try:
        rel_path = filepath.relative_to(src_root)
    except ValueError:
        # If file is not in src/, use full path
        rel_path = filepath
    
    # Convert path to components
    parts = list(rel_path.parts)
    
    # Replace filename with its stem + _H_
    file_stem = rel_path.stem
    parts[-1] = file_stem
    
    # Join with underscores and convert to uppercase
    guard_name = "_".join(parts).upper()
    
    # Replace non-alphanumeric characters with underscores
    guard_name = re.sub(r'[^A-Z0-9_]', '_', guard_name)
    
    # Add project prefix if provided
    if project_name:
        guard_name = f"{project_name.upper()}_{guard_name}"
    
    # Ensure it ends with _H_
    if not guard_name.endswith("_H_"):
        guard_name += "_H_"
    
    return guard_name


def find_include_guards(content: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Find existing include guards in the file content.
    
    Returns:
        Tuple of (ifndef_line, define_line, endif_line) or (None, None, None)
    """
    # Pattern for #ifndef - more permissive to handle various cases
    # Allows uppercase, lowercase, digits, and underscores
    # Also handles trailing comments and whitespace
    ifndef_match = re.search(r'^\s*#\s*ifndef\s+([A-Za-z0-9_]+)\s*(?://.*)?$', content, re.MULTILINE)
    
    if not ifndef_match:
        return None, None, None
    
    guard_name = ifndef_match.group(1)
    
    # Pattern for #define (should follow #ifndef, allows comments and whitespace)
    define_pattern = rf'^\s*#\s*define\s+{re.escape(guard_name)}\s*(?://.*)?$'
    define_match = re.search(define_pattern, content, re.MULTILINE)
    
    # Pattern for #endif with optional comment
    endif_pattern = rf'^\s*#\s*endif\s*(?://.*{re.escape(guard_name)}.*)?(?://.*)?$'
    endif_match = re.search(endif_pattern, content, re.MULTILINE)
    
    if ifndef_match and define_match and endif_match:
        return ifndef_match.group(0), define_match.group(0), endif_match.group(0)
    
    return None, None, None


def fix_include_guards(filepath: Path, src_root: Path, project_name: str = "", dry_run: bool = False) -> bool:
    """
    Fix the include guards of a file.
    
    Returns:
        True if the file was modified, False otherwise
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        print(f"❌ Error reading {filepath}: {e}")
        return False
    
    # Find existing guards
    ifndef_line, define_line, endif_line = find_include_guards(content)
    
    if not ifndef_line:
        print(f"⚠️  No include guard found in {filepath}")
        print(f" Add incluide guards in {filepath}")
        return add_include_guards(filepath, src_root, project_name, dry_run)
    
    # Extract current guard name
    current_guard = re.search(r'#\s*ifndef\s+([A-Za-z0-9_]+)', ifndef_line).group(1)
    
    # Generate correct name
    correct_guard = generate_guard_name(filepath, src_root, project_name)
    
    if current_guard == correct_guard:
        print(f"✅ {filepath.relative_to(src_root.parent)}: OK")
        return False
    
    print(f"🔧 {filepath.relative_to(src_root.parent)}:")
    print(f"   Old: {current_guard}")
    print(f"   New: {correct_guard}")
    
    if dry_run:
        return True
    
    # Replace the guard
    new_content = content
    
    # Replace #ifndef
    new_ifndef = f"#ifndef {correct_guard}"
    new_content = new_content.replace(ifndef_line, new_ifndef)
    
    # Replace #define
    new_define = f"#define {correct_guard}"
    new_content = new_content.replace(define_line, new_define)
    
    # Replace #endif with comment (find the last one)
    new_endif = f"#endif  // {correct_guard}"
    # Find all #endif occurrences
    endif_pattern = r'^\s*#\s*endif\s*(?://.*)?$'
    matches = list(re.finditer(endif_pattern, new_content, re.MULTILINE))
    
    if matches:
        # Replace the last occurrence
        last_match = matches[-1]
        new_content = (
            new_content[:last_match.start()] + 
            new_endif + 
            new_content[last_match.end():]
        )
    
    # Write the modified file
    try:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)
        return True
    except Exception as e:
        print(f"❌ Error writing {filepath}: {e}")
        return False

def add_include_guards(filepath: Path, src_root: Path, project_name: str = "", dry_run: bool = False) -> bool:
    """
    Add include guards to a file that doesn't have them.
    
    Returns:
        True if the file was modified, False otherwise
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        print(f"❌ Error reading {filepath}: {e}")
        return False

    guard_name = generate_guard_name(filepath, src_root, project_name)

    print(f" {filepath.relative_to(src_root.parent)}:")
    print(f"   Adding guard: {guard_name}")

    if dry_run:
        return True

    new_content = f"#ifndef {guard_name}\n#define {guard_name}\n\n{content.rstrip()}\n\n#endif  // {guard_name}\n"

    try:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)
        return True
    except Exception as e:
        print(f"❌ Error writing {filepath}: {e}")
        return False

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Fix include guards according to Google C++ style guide",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python fix_include_guards.py src/                                    # Fix all files in src/
  python fix_include_guards.py src/ --dry-run                          # Simulation mode
  python fix_include_guards.py src/ --project FUNTIDES                 # With project prefix
  python fix_include_guards.py src/ --exclude external third_party     # Exclude directories
  python fix_include_guards.py src/ --exclude external --exclude build # Exclude multiple directories
        """
    )
    
    parser.add_argument('src_dir', type=str, help='Path to the src folder')
    parser.add_argument('--project', type=str, default='', help='Project name (guard prefix)')
    parser.add_argument('--dry-run', action='store_true', help='Simulation mode (do not modify files)')
    parser.add_argument('--extensions', type=str, default='h,hpp', 
                       help='File extensions to process (comma-separated)')
    parser.add_argument('--exclude', action='append', default=[],
                       help='Directory/directories to exclude (can be used multiple times)')
    
    args = parser.parse_args()
    
    src_root = Path(args.src_dir).resolve()
    
    if not src_root.exists() or not src_root.is_dir():
        print(f"❌ Directory {src_root} does not exist")
        sys.exit(1)
    
    extensions = [f".{ext.strip()}" for ext in args.extensions.split(',')]
    exclude_dirs = args.exclude

    exclude_dirs = [d.strip("/") for d in args.exclude]
    
    print(f"🔍 Searching for header files in {src_root}")
    print(f"   Extensions: {', '.join(extensions)}")
    if args.project:
        print(f"   Project: {args.project}")
    if exclude_dirs:
        print(f"   Excluded directories: {', '.join(exclude_dirs)}")
    if args.dry_run:
        print("   Mode: DRY RUN (no modifications)")
    print()
    
    # Find all header files
    header_files = []
    for ext in extensions:
        for file in src_root.rglob(f"*{ext}"):
            if not should_exclude_path(file, src_root, exclude_dirs):
                header_files.append(file)
    
    if not header_files:
        print(f"⚠️  No header files found in {src_root}")
        sys.exit(0)
    
    print(f"📁 {len(header_files)} file(s) found\n")
    
    # Process each file
    modified_count = 0
    for filepath in sorted(header_files):
        if fix_include_guards(filepath, src_root, args.project, args.dry_run):
            modified_count += 1
    
    print(f"\n📊 Summary:")
    print(f"   Total: {len(header_files)} file(s)")
    print(f"   Modified: {modified_count} file(s)")
    
    if args.dry_run and modified_count > 0:
        print(f"\n💡 Run without --dry-run to apply changes")


if __name__ == '__main__':
    main()
