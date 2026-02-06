#!/usr/bin/env python3
"""Generate binary matrix files from test/benchmark source files."""

import struct
import re
import sys
import os

def parse_matrix_string(matrix_content):
    """Parse the matrix from the content between [ and ]."""
    # Remove string literals and concatenation
    # The format is "..." "..." with the actual content inside

    # Extract all content between quotes
    parts = re.findall(r'"([^"]*)"', matrix_content)
    full_content = ''.join(parts)

    # Remove the brackets
    full_content = full_content.strip()
    if full_content.startswith('['):
        full_content = full_content[1:]
    if full_content.endswith(']'):
        full_content = full_content[:-1]

    # Split by semicolon to get rows
    row_strs = full_content.split(';')

    rows = []
    for row_str in row_strs:
        row_str = row_str.strip()
        if not row_str:
            continue
        # Parse integers
        nums = row_str.split()
        if nums:
            row = [int(x) for x in nums]
            rows.append(row)

    return rows

def write_binary_matrix(path, rows):
    """Write matrix to binary file."""
    if not rows:
        print(f"Error: No rows to write for {path}")
        return False

    num_rows = len(rows)
    num_cols = len(rows[0])

    # Verify all rows have same length
    for i, row in enumerate(rows):
        if len(row) != num_cols:
            print(f"Error: Row {i} has {len(row)} cols, expected {num_cols}")
            return False

    print(f"Writing {path}: {num_rows} rows x {num_cols} cols = {num_rows * num_cols} values")

    with open(path, 'wb') as f:
        # Write header: rows, cols as int64
        f.write(struct.pack('<q', num_rows))
        f.write(struct.pack('<q', num_cols))

        # Write data row by row
        for row in rows:
            for val in row:
                f.write(struct.pack('<q', val))

    return True

def find_matrix_block(content, start_marker):
    """Find a matrix block after a marker, handling the tableau{ pattern."""
    pos = content.find(start_marker)
    if pos < 0:
        return None, None

    # Look for 'tableau{' or similar
    tableau_pos = content.find('tableau{', pos)
    if tableau_pos < 0 or tableau_pos > pos + 500:
        # Try 'tableau {' with space
        tableau_pos = content.find('tableau {', pos)
        if tableau_pos < 0 or tableau_pos > pos + 500:
            return None, None

    # Find the opening quote with bracket
    start = content.find('"[', tableau_pos)
    if start < 0:
        return None, None

    # Find _mat}; which ends the matrix
    end = content.find('_mat};', start)
    if end < 0:
        end = content.find('_mat}', start)
    if end < 0:
        return None, None

    # Include the _mat part
    end = end + 4  # include '_mat'

    return start, end

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_dir = os.path.dirname(script_dir)
    data_dir = os.path.join(repo_dir, 'data')

    # Create data directory if it doesn't exist
    os.makedirs(data_dir, exist_ok=True)

    # Read test file
    test_file = os.path.join(repo_dir, 'test', 'simplex_lexmin_test.cpp')
    with open(test_file, 'r') as f:
        content = f.read()

    # Extract first matrix (36x217)
    print("Extracting matrix 1 (36x217)...")
    start1, end1 = find_matrix_block(content, '"LexMinSimplexTest"_test')
    if start1 is None:
        print("Failed to find first matrix block")
        return 1

    matrix1_content = content[start1:end1]
    rows1 = parse_matrix_string(matrix1_content)
    if not rows1:
        print("Failed to parse first matrix")
        return 1

    output1 = os.path.join(data_dir, 'simplex_tableau_36x217.binmat')
    if not write_binary_matrix(output1, rows1):
        return 1

    # Extract second matrix (40x117)
    print("\nExtracting matrix 2 (40x117)...")
    start2, end2 = find_matrix_block(content, '"LexMinSimplexTest2"_test')
    if start2 is None:
        print("Failed to find second matrix block")
        return 1

    matrix2_content = content[start2:end2]
    rows2 = parse_matrix_string(matrix2_content)
    if not rows2:
        print("Failed to parse second matrix")
        return 1

    output2 = os.path.join(data_dir, 'simplex_tableau_40x117.binmat')
    if not write_binary_matrix(output2, rows2):
        return 1

    print("\nBinary matrix files generated successfully!")
    return 0

if __name__ == '__main__':
    sys.exit(main())
