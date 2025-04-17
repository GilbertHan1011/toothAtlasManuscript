def run_trajectory_conservation_analysis(
    adata_path,
    output_dir,
    pseudo_col='pseudo',
    batch_col='Sample',
    n_bins=100,
    adaptive_kernel=True,
    gene_thred=0.1,
    batch_thred=0.4,
    tail_num=0.05,
    ensure_tail=True,
    dtw_radius=3,
    use_fastdtw=True,
    normalize='zscore',
    variation_filter_level='basic',
    top_n_genes=4000,
    spline_smoothing=2,
    interpolation_factor=1,
    n_jobs=-1,
    save_figures=True,
    gene_subset=None,
    layer="logcounts"
):
    """
    Run trajectory conservation analysis on scRNA-seq data.
    
    Parameters:
    -----------
    adata_path : str
        Path to the AnnData h5ad file
    output_dir : str or Path
        Directory to save output files
    pseudo_col : str, default='pseudo'
        Column in adata.obs containing pseudotime values
    batch_col : str, default='Sample'
        Column in adata.obs containing batch/sample information
    n_bins : int, default=100
        Number of interpolation points along pseudotime
    adaptive_kernel : bool, default=True
        Whether to use adaptive kernel width for interpolation
    gene_thred : float, default=0.1
        Filter genes expressed in at least this fraction of bins
    batch_thred : float, default=0.3
        Filter batches covering at least this fraction of timeline
    ensure_tail : bool, default=True
        Ensure batches cover the tail region
    dtw_radius : int, default=3
        Radius parameter for fastdtw algorithm
    use_fastdtw : bool, default=True
        Whether to use fastdtw algorithm
    normalize : str, default='zscore'
        Method to normalize trajectories before DTW calculation
    variation_filter_level : str, default='basic'
        Level of filtering for sample variation ('off', 'basic', 'stringent')
    top_n_genes : int, default=4000
        Number of top conserved genes to select for fitting
    spline_smoothing : float, default=2
        Smoothing parameter for spline fitting
    interpolation_factor : int, default=1
        Factor for interpolation when fitting
    n_jobs : int, default=-1
        Number of parallel jobs (-1 for all available cores)
    save_figures : bool, default=True
        Whether to save visualization figures
    gene_subset : list, default=None
        Optional list of genes to use (if None, all genes are used)
    layer : str, default="logcounts"
        Layer in anndata to use for expression values
        
    Returns:
    --------
    dict
        Dictionary containing:
        - conservation_results: Results from conservation analysis
        - fit_results: Results from trajectory fitting
        - filtered_genes: List of filtered genes
        - selected_genes: List of selected genes for fitting
        - reshaped_data: 3D matrix of expression data
    """
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from pathlib import Path
    import os
    from datetime import datetime
    
    # Import from the trajDTW package
    from trajDTW import (
        anndata_to_3d_matrix, 
        calculate_trajectory_conservation,
        TrajectoryFitter,
        get_most_conserved_samples,
        fit_with_conserved_samples,
        extract_pairwise_distances,
        create_gene_position_mapping
    )
    
    # Convert output_dir to Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define variation filtering parameters
    VARIATION_FILTERING = {
        'off': {
            'filter_samples_by_variation': False
        },
        'basic': {
            'filter_samples_by_variation': True,
            'variation_threshold': 0.1,  # Minimum coefficient of variation
            'variation_metric': 'max',
            'min_valid_samples': 2       # At least 2 samples needed
        },
        'stringent': {
            'filter_samples_by_variation': True,
            'variation_threshold': 0.2, 
            'variation_metric': 'max',
            'min_valid_samples': 2
        }
    }
    
    # Setup logging
    log_file = output_dir / f"run_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    def log(message):
        """Log message to both console and file"""
        print(message)
        with open(log_file, 'a') as f:
            f.write(f"{message}\n")
    
    # Start the analysis
    log(f"\n=== Trajectory Conservation Analysis Pipeline ===")
    log(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"Input data: {adata_path}")
    log(f"Output directory: {output_dir}")
    
    # ================ 1. BUILD 3D MATRIX ================
    log("\n1. Building 3D Matrix from AnnData")
    log("-" * 50)
    
    # Load AnnData
    log("Loading AnnData...")
    adata = sc.read_h5ad(adata_path)
    
    # Apply gene subset if provided
    if gene_subset is not None:
        adata = adata[:, gene_subset].copy()
        log(f"Using subset of {len(gene_subset)} genes")
    
    log(f"AnnData shape: {adata.shape}")
    
    # Convert to 3D matrix
    log("\nConverting to 3D matrix using Gaussian kernel interpolation...")
    result = anndata_to_3d_matrix(
        adata=adata,
        pseudo_col=pseudo_col,
        batch_col=batch_col,
        n_bins=n_bins,
        adaptive_kernel=adaptive_kernel,
        gene_thred=gene_thred,
        batch_thred=batch_thred,
        tail_num= tail_num,
        ensure_tail=ensure_tail,
        layer=layer
    )
    reshaped_data = result["reshaped_data"]
    filtered_genes = result['filtered_genes']
    batch_names = result['batch_names']
    
    # Save reshaped_data as npy file
    log("Saving reshaped data as npy file...")
    np.save(output_dir / "reshaped_data.npy", reshaped_data)
    pd.DataFrame(filtered_genes).to_csv(output_dir / "filtered_genes.csv")
    pd.DataFrame(batch_names).to_csv(output_dir / "batch_names.csv")
    log(f"Reshaped data dimensions: {reshaped_data.shape}")
    log(f"Number of filtered genes: {len(filtered_genes)}")
    log(f"Batches included: {batch_names}")
    
    # ================ 2. CALCULATE CONSERVATION ================
    log("\n2. Calculating Trajectory Conservation")
    log("-" * 50)
    
    filter_params = VARIATION_FILTERING[variation_filter_level]
    log(f"Using variation filter level: {variation_filter_level}")
    
    log("Calculating trajectory conservation...")
    conservation_results = calculate_trajectory_conservation(
        trajectory_data=reshaped_data,
        gene_names=filtered_genes, 
        save_dir=output_dir,
        prefix="traj_conservation",
        dtw_radius=dtw_radius,
        use_fastdtw=use_fastdtw,
        normalize=normalize,
        **filter_params
    )
    
    # Extract and save pairwise distances
    log("Extracting pairwise distances...")
    pairwise_distances_df = extract_pairwise_distances(
        conservation_results, 
        output_csv=output_dir / "pairwise_distances.csv"
    )
    
    # Save conservation scores
    conservation_scores_df = conservation_results["conservation_scores"]
    conservation_scores_df.to_csv(output_dir / "conservation_scores.csv")
    log(f"Conservation scores saved to: {output_dir / 'conservation_scores.csv'}")
    
    # Select top conserved genes
    log(f"Selecting top {top_n_genes} conserved genes for fitting...")
    selected_genes = np.array(conservation_scores_df["gene"].head(n=top_n_genes))
    
    # Create gene position mapping
    log("Creating gene position mapping...")
    gene_mapping = create_gene_position_mapping(selected_genes, filtered_genes)
    
    # ================ 3. FIT TRAJECTORIES ================
    log("\n3. Fitting Gene Expression Trajectories")
    log("-" * 50)
    
    log("Fitting trajectories with conserved samples...")
    fit_res = fit_with_conserved_samples(
        reshaped_data=reshaped_data, 
        gene_names=selected_genes,
        gene_positions=gene_mapping,
        conserved_samples=conservation_results["conserved_samples"], 
        interpolation_factor=interpolation_factor,
        top_n_genes=None,  # Use all selected genes
        verbose=True, 
        spline_smoothing=spline_smoothing,
        n_jobs=n_jobs
    )
    
    # Save fitted trajectories
    log("Saving fitted trajectories...")
    fitdf = pd.DataFrame(fit_res["standard_results"]["fitted_trajectories"])
    fitdf.columns = fit_res["top_gene_names"]
    fitdf.to_csv(output_dir / "fitted_trajectories.csv")
    
    fitdfOptimized = pd.DataFrame(fit_res["optimized_results"]["fitted_trajectories"])
    fitdfOptimized.columns = fit_res["top_gene_names"]
    fitdfOptimized.to_csv(output_dir / "fitted_trajectories_optimized.csv")
    
    # ================ 4. VISUALIZATIONS ================
    if save_figures:
        log("\n4. Creating Visualizations")
        log("-" * 50)
        
        # 1. Heatmap of fitted trajectories
        log("Creating heatmap of fitted trajectories...")
        plt.figure(figsize=(16, 12))
        sns.heatmap(fitdf, cmap='viridis', cbar=True)
        plt.title("Fitted Gene Expression Trajectories")
        plt.xlabel("Genes")
        plt.ylabel("Pseudotime Points")
        plt.savefig(output_dir / "fitted_trajectories_heatmap.png", bbox_inches='tight', dpi=300)
        
        # 2. K-means clustering of genes
        log("Performing K-means clustering of genes...")
        from sklearn.cluster import KMeans
        
        # Choose number of clusters (this can be parameterized)
        n_clusters = 8
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(fitdf.T)  # Transpose to cluster genes
        
        # Add cluster labels to gene data
        gene_clusters = pd.DataFrame({
            'gene': fit_res["top_gene_names"],
            'cluster': cluster_labels
        })
        gene_clusters.to_csv(output_dir / "gene_clusters.csv", index=False)
        
        # 3. Clustered heatmap
        log("Creating clustered heatmap...")
        # Sort genes by cluster
        order = gene_clusters.sort_values('cluster').index
        fitdf_clustered = fitdf[fit_res["top_gene_names"][order]]
        
        plt.figure(figsize=(16, 12))
        sns.heatmap(fitdf_clustered, cmap='viridis', cbar=True)
        plt.title("K-means Clustered Gene Expression Trajectories")
        plt.xlabel("Genes (Clustered)")
        plt.ylabel("Pseudotime Points")
        plt.savefig(output_dir / "clustered_fitted_trajectories_heatmap.png", bbox_inches='tight', dpi=300)
        
        # 4. Cluster profile plot
        log("Creating cluster profile plot...")
        plt.figure(figsize=(14, 10))
        
        for cluster in range(n_clusters):
            cluster_genes = gene_clusters[gene_clusters['cluster'] == cluster]['gene']
            if len(cluster_genes) > 0:
                cluster_data = fitdf[cluster_genes].mean(axis=1)
                plt.plot(cluster_data, label=f'Cluster {cluster} (n={len(cluster_genes)})')
        
        plt.title("Average Expression Profiles by Cluster")
        plt.xlabel("Pseudotime Points")
        plt.ylabel("Average Expression")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(output_dir / "cluster_profiles.png", bbox_inches='tight', dpi=300)
    
    log("\n=== Analysis Complete ===")
    log(f"Finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Return results
    return {
        'conservation_results': conservation_results,
        'fit_results': fit_res,
        'filtered_genes': filtered_genes,
        'selected_genes': selected_genes,
        'reshaped_data': reshaped_data
    }