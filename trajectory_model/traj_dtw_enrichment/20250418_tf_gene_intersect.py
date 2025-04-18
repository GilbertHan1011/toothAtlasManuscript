#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import sys
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

def read_tf_gene_data(file_path):
    """
    Read TF-gene data from CSV file into pandas DataFrame.
    
    Args:
        file_path (str): Path to CSV file containing TF, Gene, Score columns
        
    Returns:
        pandas.DataFrame: DataFrame with TF, Gene, Score columns
    """
    try:
        # Check if file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")
        
        # Read CSV file
        df = pd.read_csv(file_path)
        
        # Check column names
        required_cols = ['TF', 'Gene', 'Score']
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Required column '{col}' not found in {file_path}")
        
        # Ensure Score column is numeric
        df['Score'] = pd.to_numeric(df['Score'], errors='coerce')
        
        # Drop rows with NaN scores
        df = df.dropna(subset=['Score'])
        
        return df
    
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        raise

def filter_by_score(df, threshold):
    """
    Filter DataFrame to keep only rows with Score above threshold.
    
    Args:
        df (pandas.DataFrame): DataFrame with TF, Gene, Score columns
        threshold (float): Minimum score threshold
        
    Returns:
        pandas.DataFrame: Filtered DataFrame
    """
    filtered_df = df[df['Score'] >= threshold].copy()
    return filtered_df

def find_intersection(df1, df2):
    """
    Find TF-gene pairs that are common between two DataFrames.
    
    Args:
        df1 (pandas.DataFrame): First DataFrame with TF, Gene, Score columns
        df2 (pandas.DataFrame): Second DataFrame with TF, Gene, Score columns
        
    Returns:
        pandas.DataFrame: DataFrame with common TF-Gene pairs and scores from both files
    """
    # Create composite key for TF-Gene pair
    df1['TF_Gene'] = df1['TF'] + '|' + df1['Gene']
    df2['TF_Gene'] = df2['TF'] + '|' + df2['Gene']
    
    # Find common TF-Gene pairs
    common_pairs = set(df1['TF_Gene']).intersection(set(df2['TF_Gene']))
    
    if not common_pairs:
        print("No common TF-Gene pairs found between the two files.")
        return pd.DataFrame(columns=['TF', 'Gene'])
    
    # Filter to keep only common pairs
    df1_common = df1[df1['TF_Gene'].isin(common_pairs)].copy()
    
    # Only keep TF and Gene columns, drop duplicates
    result = df1_common[['TF', 'Gene']].drop_duplicates().reset_index(drop=True)
    
    return result

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Find intersection of TF-gene pairs from two files with score filtering.')
    parser.add_argument('--file1', '-f1', required=True,
                       help='Path to first CSV file containing TF, Gene, Score columns')
    parser.add_argument('--file2', '-f2', required=True,
                       help='Path to second CSV file containing TF, Gene, Score columns')
    parser.add_argument('--threshold1', '-t1', type=float, default=0.0,
                       help='Score threshold for first file (default: 0.0)')
    parser.add_argument('--threshold2', '-t2', type=float, default=0.0,
                       help='Score threshold for second file (default: 0.0)')
    parser.add_argument('--output', '-o', required=True,
                       help='Path to output CSV file for intersection results')
    
    args = parser.parse_args()
    
    try:
        # Read input files
        print(f"Reading first file: {args.file1}")
        df1 = read_tf_gene_data(args.file1)
        print(f"Found {len(df1)} TF-gene pairs in first file")
        
        print(f"Reading second file: {args.file2}")
        df2 = read_tf_gene_data(args.file2)
        print(f"Found {len(df2)} TF-gene pairs in second file")
        
        # Filter by score thresholds
        print(f"Filtering first file with threshold {args.threshold1}")
        df1_filtered = filter_by_score(df1, args.threshold1)
        print(f"Kept {len(df1_filtered)} TF-gene pairs after filtering")
        
        print(f"Filtering second file with threshold {args.threshold2}")
        df2_filtered = filter_by_score(df2, args.threshold2)
        print(f"Kept {len(df2_filtered)} TF-gene pairs after filtering")
        
        # Find intersection
        print("Finding intersection of TF-gene pairs")
        intersection_df = find_intersection(df1_filtered, df2_filtered)
        print(f"Found {len(intersection_df)} common TF-gene pairs")
        
        # Ensure output directory exists
        output_path = expand_path(args.output)
        output_dir = os.path.dirname(output_path)
        
        if output_dir and not os.path.exists(output_dir):
            print(f"Creating output directory: {output_dir}")
            os.makedirs(output_dir, exist_ok=True)
        
        # Write intersection to output file
        intersection_df.to_csv(output_path, index=False)
        print(f"Wrote intersection to: {output_path}")
        
        # Print sample of intersection
        print("\nSample of intersection (first 5 rows):")
        print(intersection_df.head(5).to_string(index=False))
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 