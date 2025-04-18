#!/usr/bin/env python3
import csv
import ast
import os
import sys
import argparse
import re
from pathlib import Path

def expand_path(path):
    """
    Expand user home directory symbol (~) and return an absolute path.
    
    Args:
        path (str): Path that may contain ~ or relative components
        
    Returns:
        str: Absolute path with ~ expanded
    """
    return os.path.abspath(os.path.expanduser(path))

def is_list_of_tuples(s):
    """
    Check if a string appears to represent a list of tuples with gene names and scores.
    
    Args:
        s (str): String to check
        
    Returns:
        bool: True if the string appears to represent a list of tuples
    """
    return s.strip().startswith('[') and s.strip().endswith(']') and '(' in s and ')' in s

def parse_target_genes(target_genes_str):
    """
    Parse target genes string into a list of (gene, score) tuples.
    Handles different formats that might appear in the data.
    
    Args:
        target_genes_str (str): String representation of target genes
        
    Returns:
        list: List of (gene, score) tuples
    """
    # Check if this appears to be a list of tuples
    if is_list_of_tuples(target_genes_str):
        try:
            # Standard format - list of tuples
            return ast.literal_eval(target_genes_str)
        except (SyntaxError, ValueError) as e:
            # If ast.literal_eval fails, we'll try regex parsing as a fallback
            pass
    
    # If we're here, either the format wasn't a list of tuples or ast.literal_eval failed
    # Try to extract gene names and scores using regex
    results = []
    
    # Pattern for matching ('GeneName', score) tuples
    pattern = r"\('([^']+)',\s*([\d\.]+)\)"
    matches = re.findall(pattern, target_genes_str)
    
    if matches:
        for gene, score in matches:
            try:
                results.append((gene, float(score)))
            except ValueError:
                # Skip if score can't be converted to float
                continue
    
    return results

def parse_bone_regulon(file_path, max_rows=None):
    """
    Parse the bone regulon CSV file to extract TF, gene and score relationships.
    
    Args:
        file_path (str): Path to the regulon CSV file
        max_rows (int, optional): Maximum number of rows to process, for testing
        
    Returns:
        list: List of tuples containing (tf, gene, score) relationships
    """
    print(f"Parsing regulon file: {file_path}")
    if max_rows:
        print(f"Limiting to {max_rows} rows for testing")
        
    file_path = expand_path(file_path)
    
    # Check if file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    # Initialize result list
    results = []
    error_count = 0
    row_count = 0
    processed_count = 0
    
    try:
        # Open the CSV file
        with open(file_path, 'r') as f:
            reader = csv.reader(f)
            
            # Skip header rows (first 3 rows)
            for _ in range(3):
                next(reader)
            
            # Process each row
            for row in reader:
                row_count += 1
                
                # Check if we've reached the maximum rows
                if max_rows and processed_count >= max_rows:
                    print(f"Reached maximum row limit of {max_rows}")
                    break
                    
                if not row:  # Skip empty rows
                    continue
                
                processed_count += 1
                tf = row[0]  # Transcription factor
                
                # Only process if we have a TF (skip empty rows)
                if tf:
                    # Parse the target genes - it's a string representation of a list of tuples
                    # The target genes are in column 9 (index 8)
                    if len(row) > 8 and row[8]:
                        try:
                            # Parse the target genes using our custom parser
                            target_genes = parse_target_genes(row[8])
                            
                            # Each target is a tuple of (gene_name, score)
                            if target_genes:
                                for gene, score in target_genes:
                                    results.append((tf, gene, score))
                            else:
                                # If no target genes were extracted, log this
                                error_count += 1
                                print(f"Warning: No target genes could be extracted for TF {tf} in row {row_count}")
                                if isinstance(row[8], str) and len(row[8]) > 50:
                                    print(f"Target genes string starts with: {row[8][:50]}...")
                        except Exception as e:
                            error_count += 1
                            print(f"Error parsing row {row_count} for TF {tf}: {e}")
                            if isinstance(row[8], str) and len(row[8]) > 50:
                                print(f"Target genes string starts with: {row[8][:50]}...")
    
    except Exception as e:
        print(f"Unexpected error while processing file: {e}")
        raise
    
    print(f"Processed {processed_count} rows with {error_count} parsing errors")
    print(f"Extracted {len(results)} gene relationships")
    
    return results

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Extract TF, gene, and score information from bone regulon CSV.')
    parser.add_argument('--input', '-i', 
                        default='/home/gilberthan/Desktop/disk2/202409_tooth/process/trajectory/20250417_mineralization_downstream/data/bone_regulon.csv',
                        help='Path to input bone regulon CSV file')
    parser.add_argument('--output', '-o', 
                        default='/home/gilberthan/Desktop/disk2/202409_tooth/process/trajectory/20250417_mineralization_downstream/data/bone_regulon_simplified.csv',
                        help='Path to output simplified CSV file')
    parser.add_argument('--max-rows', '-m', type=int, default=None,
                        help='Maximum number of rows to process (for testing)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print detailed progress information')
    
    args = parser.parse_args()
    
    # Parse the data
    try:
        results = parse_bone_regulon(args.input, args.max_rows)
    except Exception as e:
        print(f"Error parsing regulon file: {e}")
        sys.exit(1)
    
    # If no results were found, exit
    if not results:
        print("No gene relationships were extracted. Check the file format.")
        sys.exit(1)
    
    # Ensure output directory exists
    output_path = expand_path(args.output)
    output_dir = os.path.dirname(output_path)
    
    if not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)
    
    # Print the data to stdout in CSV format
    print("\nSample of extracted data (first 5 rows):")
    print("TF,Gene,Score")
    for tf, gene, score in results[:5]:
        print(f"{tf},{gene},{score}")
    
    # Write to a new CSV file
    try:
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["TF", "Gene", "Score"])
            writer.writerows(results)
            
        print(f"\nProcessed {len(results)} gene relationships and saved to '{output_path}'")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()