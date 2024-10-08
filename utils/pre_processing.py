import numpy as np
import scanpy as sc
import anndata
import pandas as pd
import matplotlib.pyplot as plt
import collections
import glob


def add_up_duplicate_gene_name_columns(adata, print_gene_names=True, verbose=False):
    """ This function finds duplicate gene names in adata.var (i.e. duplicates 
    in the index of adata.var). For each cell, it adds up the counts of columns 
    with the same gene name, removes the duplicate columns, and replaces the 
    counts of the remaining column with the sum of the duplicates.
    Returns anndata object."""
    print("TO DO: STORE ENSEMBL IDS OF MERGED IDS AND MATCHING COUNTS IN A SEPARATE COLUMN!!!")
    duplicate_boolean = adata.var.index.duplicated()
    duplicate_genes = adata.var.index[duplicate_boolean]
    print("Number of duplicate genes: " + str(len(duplicate_genes)))
    if print_gene_names == True:
        print(duplicate_genes)
    columns_to_replace = list()
    columns_to_remove = list()
    new_columns_array = np.empty((adata.shape[0], 0))
    for gene in duplicate_genes:
        if verbose:
            print("Calculating for gene", gene)
        # get rows in adata.var with indexname equal to gene
        # indexing zero here is to get rid of tuple output and access array
        gene_colnumbers = np.where(adata.var.index == gene)[0]
        # isolate the columns
        gene_counts = adata.X[:, gene_colnumbers].copy()
        # add up gene counts and add new column to new_columns_array
        new_columns_array = np.hstack((new_columns_array, np.sum(gene_counts, axis=1)))
        # store matching column location in real adata in list:
        columns_to_replace.append(gene_colnumbers[0])
        # store remaining column locations in columns to remove:
        columns_to_remove = columns_to_remove + gene_colnumbers[1:].tolist()
    # replace first gene column with new col:
    adata.X[:, columns_to_replace] = new_columns_array
    # remove remaining duplicate columns:
    columns_to_keep = [
        i for i in np.arange(adata.shape[1]) if i not in columns_to_remove
    ]
    adata = adata[:, columns_to_keep].copy()
    if verbose:
        print("Done!")
    return adata
    

def add_cell_annotations(adata, var_index):
    """ This function adds annotation to anndata:  
    cell level:  
    total_counts, log10_total_counts, n_genes_detected, mito_frac, ribo_frac,   
    compl(exity)  
    gene_level:  
    n_cells 

    Arguments:
        adata - anndata object, raw (unnormalized!)
        var_index < "gene_symbols", "gene_ids" > - set to which type of gene
            naming is used in adata.var.index

    Returns:
        anndata object with annotations  
    """
    # cells:
    # total transcript count per cell
    adata.obs['total_counts'] = np.sum(adata.X, axis=1)
    adata.obs['log10_total_counts'] = np.log10(adata.obs['total_counts'])
    # number of genes expressed
    # translate matrix to boolean (True if count is larger than 0):
    boolean_expression = adata.X > 0
    adata.obs['n_genes_detected'] = np.sum(boolean_expression, axis=1)
    # fraction mitochondrial transcripts per cell
    if var_index == "gene_symbols":
        mito_genes = [
            gene for gene in adata.var.index if gene.lower().startswith("mt-")
        ]
    elif var_index == "gene_ids":
        mito_genes = [
            gene_id
            for gene_id, gene_symbol in zip(adata.var.index, adata.var.gene_symbols)
            if gene_symbol.lower().startswith("mt-")
        ]
    else:
        raise ValueError(
            "var_index argument should be set to either gene_symbols or gene_ids"
        )
    # conversion to array in line below is necessary if adata.X is sparse
    adata.obs['mito_frac'] = np.array(np.sum(
        adata[:,mito_genes].X, axis=1)).flatten() / adata.obs['total_counts']
    # fraction ribosomal transcripts per cell
    if var_index == "gene_symbols":
        ribo_genes = [
            gene
            for gene in adata.var.index
            if (gene.lower().startswith("rpl") or gene.lower().startswith("rps"))
        ]
    elif var_index == "gene_ids":
        ribo_genes = [
            gene_id
            for gene_id, gene_symbol in zip(adata.var.index, adata.var.gene_symbols)
            if (
                gene_symbol.lower().startswith("rpl")
                or gene_symbol.lower().startswith("rps")
            )
        ]
    adata.obs['ribo_frac'] = np.array(np.sum(
        adata[:,ribo_genes].X, axis=1)).flatten() / adata.obs['total_counts']
    # cell complexity (i.e. number of genes detected / total transcripts)
    adata.obs['compl'] = adata.obs['n_genes_detected']\
    / adata.obs['total_counts']
    # genes
    adata.var['n_cells'] = np.sum(boolean_expression, axis=0).T
    return adata

def merge_metadata(adata, meta):
    """
    Merge adata.obs with another DataFrame (meta) based on index.
    
    Parameters:
    adata: An object with an 'obs' attribute that is a DataFrame.
    meta: A DataFrame to merge with adata.obs.
    
    Returns:
    None: Updates adata.obs in place.
    """
    A = adata.obs
    B = meta
    
    # Check for overlapping column names
    overlapping_cols = A.columns.intersection(B.columns)

    if not overlapping_cols.empty:
        print(f"Warning: Columns {overlapping_cols.tolist()} exist in both DataFrames. They will be replaced.")
        # Merge and replace overlapping columns
        result = A.combine_first(B)
    else:
        print("No overlapping columns found. Merging without replacement.")
        # Merge without replacement
        result = A.join(B)

    # Update adata.obs with the merged result
    adata.obs = result
    return(adata)
