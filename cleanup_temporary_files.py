#!/usr/bin/env python3
"""
Standalone script to clean up temporary pkl files from PACE analysis directories.
This can be used to clean up existing directories that have accumulated temporary files.

Usage:
    python cleanup_temporary_files.py /path/to/output/directory
    python cleanup_temporary_files.py /path/to/output/directory --dry-run
"""

import os
import sys
import glob
import argparse


def cleanup_temporary_pkl_files(output_dir, dry_run=False):
    """
    Remove temporary pkl files to save disk space.
    
    Args:
        output_dir: Directory to clean up
        dry_run: If True, only show what would be deleted without actually deleting
    """
    if not os.path.exists(output_dir):
        print(f"Error: Directory {output_dir} does not exist")
        return
    
    files_to_remove = []
    total_size = 0
    
    # Define patterns for temporary files to remove
    # Use os.listdir to avoid glob pattern issues with special characters in filenames
    try:
        all_files = os.listdir(output_dir)
        for filename in all_files:
            file_path = os.path.join(output_dir, filename)
            if os.path.isfile(file_path):
                # Check if this file matches our cleanup criteria
                should_remove = False
                
                if (filename.endswith('_threshold_dict.pkl') or      # Threshold dictionary files
                    filename.startswith('di_') and filename.endswith('.pkl') or  # Di dictionary files
                    filename.startswith('performance_dict_') and filename.endswith('.pkl') or  # Performance dictionary files
                    filename.startswith('all_par_') and filename.endswith('.pkl')):  # All parameters files
                    should_remove = True
                
                if should_remove:
                    file_size = os.path.getsize(file_path)
                    files_to_remove.append((file_path, file_size))
                    total_size += file_size
    except Exception as e:
        print(f"Error reading directory {output_dir}: {e}")
        return
    
    # Sort by size (largest first) for better reporting
    files_to_remove.sort(key=lambda x: x[1], reverse=True)
    
    if not files_to_remove:
        print(f"No temporary pkl files found in {output_dir}")
        return
    
    # Show summary
    size_mb = total_size / (1024 * 1024)
    print(f"Found {len(files_to_remove)} temporary pkl files totaling {size_mb:.1f} MB")
    
    if dry_run:
        print("\nDRY RUN - Files that would be removed:")
        for file_path, file_size in files_to_remove:
            size_mb = file_size / (1024 * 1024)
            print(f"  {os.path.basename(file_path):50} ({size_mb:6.1f} MB)")
        print(f"\nTotal: {len(files_to_remove)} files, {total_size / (1024 * 1024):.1f} MB")
        print("Run without --dry-run to actually delete these files.")
        return
    
    # Remove files
    removed_count = 0
    removed_size = 0
    errors = []
    
    print("\nRemoving files:")
    for file_path, file_size in files_to_remove:
        try:
            os.remove(file_path)
            removed_count += 1
            removed_size += file_size
            size_mb = file_size / (1024 * 1024)
            print(f"  ✓ {os.path.basename(file_path):50} ({size_mb:6.1f} MB)")
        except Exception as e:
            errors.append(f"  ✗ {os.path.basename(file_path)}: {e}")
    
    # Summary
    removed_mb = removed_size / (1024 * 1024)
    print(f"\nSuccessfully removed {removed_count} files, freed {removed_mb:.1f} MB")
    
    if errors:
        print(f"\nErrors encountered ({len(errors)} files):")
        for error in errors:
            print(error)


def main():
    parser = argparse.ArgumentParser(
        description="Clean up temporary pkl files from PACE analysis directories"
    )
    parser.add_argument(
        "directory",
        help="Path to the output directory to clean up"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be deleted without actually deleting"
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Clean up all subdirectories recursively"
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.directory):
        print(f"Error: Directory {args.directory} does not exist")
        sys.exit(1)
    
    if args.recursive:
        # Find all subdirectories and clean them up
        total_removed = 0
        total_size = 0
        
        for root, dirs, files in os.walk(args.directory):
            # Skip hidden directories and common non-analysis directories
            dirs[:] = [d for d in dirs if not d.startswith('.') and d not in ['__pycache__']]
            
            # Check if this directory has any pkl files
            pkl_files = glob.glob(os.path.join(root, "*.pkl"))
            if pkl_files:
                print(f"\nCleaning directory: {root}")
                cleanup_temporary_pkl_files(root, args.dry_run)
    else:
        cleanup_temporary_pkl_files(args.directory, args.dry_run)
        
        # Also check immediate subdirectories (common case for bias correction runs)
        if not args.dry_run:
            print(f"\nChecking subdirectories in {args.directory}...")
        
        try:
            for item in os.listdir(args.directory):
                item_path = os.path.join(args.directory, item)
                if os.path.isdir(item_path) and not item.startswith('.'):
                    # Check if this subdirectory has pkl files
                    pkl_files = glob.glob(os.path.join(item_path, "*.pkl"))
                    if pkl_files:
                        if not args.dry_run:
                            print(f"\nCleaning subdirectory: {item_path}")
                        cleanup_temporary_pkl_files(item_path, args.dry_run)
        except Exception as e:
            print(f"Error checking subdirectories: {e}")


if __name__ == "__main__":
    main() 