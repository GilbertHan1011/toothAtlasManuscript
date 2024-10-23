import pandas as pd
import argparse
import sys

# This Python script updates and merges metadata from two CSV files: a general metadata file and a single-cell metadata file. Using the pandas library, it cleans the data by removing unnecessary columns and handling missing values. The script performs a left join on the "Sample" column to combine the datasets and saves the updated metadata to a new CSV file.
def update_metadata(meta_path, scMeta_path, output_path):
    # Read the metadata files
    meta = pd.read_csv(meta_path, index_col=0)
    scMeta = pd.read_csv(scMeta_path, index_col=0)

    # Drop columns that start with 'Unnamed'
    meta = meta.drop(columns=meta.filter(like='Unnamed').columns)
    
    # Drop the index with NaN
    meta = meta[~meta.index.isna()]
    
    # Temporarily replace NaN values
    meta = meta.fillna("Not provided")
    
    
    # Reset index of scMeta for merging
    scMeta = scMeta.reset_index()
    
    # Perform a left join on the "Sample" column
    metaUpdate = pd.merge(scMeta, meta, on="Sample", how="left")
    
    # Set the original index back
    metaUpdate.set_index('index', inplace=True)
    print(f"writing to {output_path}")
    # Save the updated DataFrame to a CSV file
    metaUpdate.to_csv(output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update metadata based on provided paths.")
    
    parser.add_argument('-f1', '--meta_path', type=str, 
                        default="../../data/metadata/metadata_latest.csv", 
                        help="Path to the metadata file (default: ../../processed_data/metadata/scMetadata_latest.csv)")
    
    parser.add_argument('-f2', '--scMeta_path', type=str, 
                        default="../../processed_data/metadata/scMetadata_latest.csv", 
                        help="Path to the single-cell metadata file.")
    
    parser.add_argument('-o', '--output_path', type=str,
                        default="../../processed_data/metadata/scMetadata_latest.csv", 
                        help="Path to save the updated metadata file.")
    
    args = parser.parse_args()

    # Call the update_metadata function with the provided arguments
    update_metadata(args.meta_path, args.scMeta_path, args.output_path)
