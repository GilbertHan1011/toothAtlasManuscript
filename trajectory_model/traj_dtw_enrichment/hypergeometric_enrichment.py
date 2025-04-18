#!/usr/bin/env python3
import argparse
import pandas as pd
from scipy.stats import hypergeom
import os

def read_gene_list(file_path):
    """Read a gene list from a file (one gene per line, or comma/space separated)."""
    genes = set()
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Split by comma, space, or tab
            for token in line.replace(',', ' ').replace('\t', ' ').split():
                genes.add(token.strip())
    return genes

def read_gmt(file_path):
    """Read a GMT file (category\tdescription\tgene1\tgene2...) and return dict of category -> set(genes)."""
    categories = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 3:
                continue
            cat = parts[0]
            genes = set(parts[2:])
            categories[cat] = genes
    return categories

def hypergeometric_enrichment(input_genes, categories, background_size=20000):
    """
    Perform hypergeometric enrichment for each category.
    Returns a list of dicts with results.
    """
    results = []
    input_genes = set(input_genes)
    N = len(input_genes)
    M = background_size
    for cat, cat_genes in categories.items():
        n = len(cat_genes)
        k = len(input_genes & cat_genes)
        # P-value: probability of >=k overlap
        pval = hypergeom.sf(k-1, M, n, N) if k > 0 else 1.0
        results.append({
            'Category': cat,
            'Category_Size': n,
            'Input_Size': N,
            'Overlap': k,
            'P-value': pval,
            'Overlapping_Genes': ','.join(sorted(input_genes & cat_genes))
        })
    return results

def main():
    parser = argparse.ArgumentParser(description='Hypergeometric enrichment analysis for a gene set.')
    parser.add_argument('--geneset', '-g', required=True, help='Input gene set file (one gene per line or comma/space separated)')
    parser.add_argument('--categories', '-c', required=True, help='Category gene sets in GMT format (category\tdescription\tgene1\tgene2...)')
    parser.add_argument('--output', '-o', required=True, help='Output CSV file for enrichment results')
    parser.add_argument('--background', '-b', type=int, default=20000, help='Background gene number (default: 20000)')
    args = parser.parse_args()

    input_genes = read_gene_list(args.geneset)
    categories = read_gmt(args.categories)
    results = hypergeometric_enrichment(input_genes, categories, background_size=args.background)
    df = pd.DataFrame(results)
    df = df.sort_values('P-value')
    df.to_csv(args.output, index=False)
    print(f"Wrote enrichment results to {args.output}")
    print(df.head(10).to_string(index=False))

if __name__ == '__main__':
    main() 