"""
Trajectory visualization module
==============================

This module provides functions for visualizing trajectory data, conservation scores,
and model fitting results.

Functions:
- plot_gene_trajectories: Plot gene expression trajectories
- plot_conservation_scores: Plot conservation scores for genes
- plot_pairwise_distances: Plot pairwise distances between trajectories
- plot_fitted_trajectories: Plot fitted trajectory models
- plot_fitting_comparison: Plot comparison of different fitting methods
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Union, Optional, Any
import warnings
import logging
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from scipy.stats import zscore

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def plot_gene_trajectories(
    data_3d: np.ndarray,
    gene_idx: int,
    gene_name: Optional[str] = None,
    sample_indices: Optional[List[int]] = None,
    pseudotime: Optional[np.ndarray] = None,
    figsize: Tuple[int, int] = (10, 6),
    color_map: str = 'viridis',
    title: Optional[str] = None,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    normalize: bool = False,
    alpha: float = 0.7,
    linewidth: float = 2.0,
    show_legend: bool = True,
    save_path: Optional[str] = None,
    ax: Optional[plt.Axes] = None
) -> plt.Figure:
    """
    Plot gene expression trajectories for a specific gene.
    
    Parameters
    ----------
    data_3d : numpy.ndarray
        3D data matrix [batches, timepoints, genes]
    gene_idx : int
        Index of the gene to plot
    gene_name : str, optional
        Name of the gene for plot title
    sample_indices : list, optional
        Indices of samples to include
    pseudotime : numpy.ndarray, optional
        Pseudotime values for x-axis (if None, uses indices)
    figsize : tuple, optional
        Figure size (width, height)
    color_map : str, optional
        Matplotlib colormap name
    title : str, optional
        Plot title (if None, uses gene_name)
    xlim : tuple, optional
        x-axis limits (min, max)
    ylim : tuple, optional
        y-axis limits (min, max)
    normalize : bool, optional
        Whether to normalize gene expression (z-score)
    alpha : float, optional
        Transparency of trajectories
    linewidth : float, optional
        Line width of trajectories
    show_legend : bool, optional
        Whether to show legend
    save_path : str, optional
        Path to save figure
    ax : matplotlib.axes.Axes, optional
        Axes to plot on (if None, creates new figure)
        
    Returns
    -------
    matplotlib.figure.Figure
        The figure with plotted trajectories
    """
    # Extract gene data
    gene_data = data_3d[:, :, gene_idx]
    
    # Filter by samples if provided
    if sample_indices is not None:
        # Validate sample indices
        valid_indices = [i for i in sample_indices if i < gene_data.shape[0]]
        if len(valid_indices) < len(sample_indices):
            warnings.warn(f"Removed {len(sample_indices) - len(valid_indices)} invalid sample indices")
        gene_data = gene_data[valid_indices]
    
    # Handle pseudotime
    x_values = np.arange(gene_data.shape[1]) if pseudotime is None else pseudotime
    
    # Normalize if requested
    if normalize:
        for i in range(gene_data.shape[0]):
            if not np.all(np.isnan(gene_data[i])):
                gene_data[i] = zscore(gene_data[i], nan_policy='omit')
    
    # Create figure if no axes provided
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    
    # Set up colormap
    cmap = cm.get_cmap(color_map)
    colors = [cmap(i / max(1, gene_data.shape[0] - 1)) for i in range(gene_data.shape[0])]
    
    # Plot each trajectory
    for i in range(gene_data.shape[0]):
        label = f"Sample {i}" if show_legend else None
        ax.plot(x_values, gene_data[i], color=colors[i], alpha=alpha, 
                linewidth=linewidth, label=label)
    
    # Set title
    if title is None:
        title = f"Expression Trajectories for {gene_name}" if gene_name else f"Gene {gene_idx}"
    ax.set_title(title)
    
    # Set labels
    ax.set_xlabel("Pseudotime" if pseudotime is not None else "Time")
    ax.set_ylabel("Expression" + (" (z-score)" if normalize else ""))
    
    # Set limits if provided
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    
    # Add legend if requested
    if show_legend:
        # Adjust for too many trajectories
        if gene_data.shape[0] > 10:
            # Use a subset for the legend
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[:10], labels[:10], loc='best', title="Samples")
        else:
            ax.legend(loc='best', title="Samples")
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Tight layout
    plt.tight_layout()
    
    # Save if path provided
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig

def plot_conservation_scores(
    conservation_results: Dict[str, Any],
    gene_names: List[str],
    top_n: int = 20,
    figsize: Tuple[int, int] = (12, 8),
    colors: Tuple[str, str] = ('#3498db', '#e74c3c'),
    show_normalized: bool = True,
    title: Optional[str] = None,
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Plot conservation scores for top genes.
    
    Parameters
    ----------
    conservation_results : dict
        Results from calculate_trajectory_conservation
    gene_names : list
        List of gene names
    top_n : int, optional
        Number of top genes to display
    figsize : tuple, optional
        Figure size (width, height)
    colors : tuple, optional
        Colors for (raw scores, normalized scores)
    show_normalized : bool, optional
        Whether to show normalized scores
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    matplotlib.figure.Figure
        The figure with plotted conservation scores
    """
    # Extract raw scores
    conservation_scores = conservation_results.get('conservation_scores')
    if conservation_scores is None:
        raise ValueError("Conservation results missing required data")
        
    # Get normalized scores if needed
    if show_normalized:
        normalized_scores = conservation_results.get('normalized_scores')
        if normalized_scores is None:
            show_normalized = False
            warnings.warn("Normalized scores not found in conservation results. Showing only raw scores.")
    
    # Check gene names
    if len(gene_names) != len(conservation_scores):
        raise ValueError(f"Length of gene_names ({len(gene_names)}) doesn't match conservation_scores ({len(conservation_scores)})")
    
    # Create dataframe
    data = {'gene': gene_names, 'score': conservation_scores}
    if show_normalized:
        data['normalized'] = normalized_scores
    df = pd.DataFrame(data)
    
    # Sort by raw score (higher is better for conservation)
    df = df.sort_values('score', ascending=False)
    
    # Get top N genes
    df_top = df.head(top_n)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot bars
    x = np.arange(len(df_top))
    bar_width = 0.35
    
    # Plot raw scores
    ax.bar(x - bar_width/2 if show_normalized else x, 
           df_top['score'], 
           width=bar_width, 
           color=colors[0], 
           label='Raw Conservation Score')
    
    # Plot normalized scores if requested
    if show_normalized:
        ax2 = ax.twinx()
        ax2.bar(x + bar_width/2, 
                df_top['normalized'], 
                width=bar_width, 
                color=colors[1], 
                alpha=0.7,
                label='Normalized Score')
        ax2.set_ylabel('Normalized Score')
        ax2.set_ylim(0, 1.05)
    
    # Set x-axis labels
    ax.set_xticks(x)
    ax.set_xticklabels(df_top['gene'], rotation=45, ha='right')
    
    # Set y-axis label
    ax.set_ylabel('Conservation Score')
    
    # Set title
    if title is None:
        title = f"Top {top_n} Conserved Genes"
    ax.set_title(title)
    
    # Add legend
    if show_normalized:
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
    else:
        ax.legend(loc='upper right')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save if path provided
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig

def plot_pairwise_distances(
    distance_matrix: np.ndarray,
    gene_name: Optional[str] = None,
    gene_idx: Optional[int] = None,
    sample_labels: Optional[List[str]] = None,
    figsize: Tuple[int, int] = (10, 8),
    cmap: str = 'viridis',
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    title: Optional[str] = None,
    colorbar_label: str = 'Distance',
    annotate: bool = True,
    annotation_fmt: str = '.2f',
    annotation_fontsize: int = 8,
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Plot pairwise distance matrix for a specific gene.
    
    Parameters
    ----------
    distance_matrix : numpy.ndarray
        Square matrix of pairwise distances
    gene_name : str, optional
        Name of the gene
    gene_idx : int, optional
        Index of the gene
    sample_labels : list, optional
        Labels for samples
    figsize : tuple, optional
        Figure size (width, height)
    cmap : str, optional
        Colormap name
    vmin : float, optional
        Minimum value for colormap
    vmax : float, optional
        Maximum value for colormap
    title : str, optional
        Plot title
    colorbar_label : str, optional
        Label for colorbar
    annotate : bool, optional
        Whether to annotate heatmap cells with values
    annotation_fmt : str, optional
        Format string for annotations
    annotation_fontsize : int, optional
        Font size for annotations
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    matplotlib.figure.Figure
        The figure with plotted pairwise distances
    """
    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("seaborn is required for plotting distance matrices")
    
    # Create sample labels if not provided
    if sample_labels is None:
        sample_labels = [f"Sample {i}" for i in range(distance_matrix.shape[0])]
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot heatmap
    sns.heatmap(
        distance_matrix,
        annot=annotate,
        fmt=annotation_fmt,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        square=True,
        xticklabels=sample_labels,
        yticklabels=sample_labels,
        cbar_kws={'label': colorbar_label},
        annot_kws={'fontsize': annotation_fontsize},
        ax=ax
    )
    
    # Set title
    if title is None:
        if gene_name:
            title = f"Pairwise Distances for {gene_name}"
        elif gene_idx is not None:
            title = f"Pairwise Distances for Gene {gene_idx}"
        else:
            title = "Pairwise Distances"
    ax.set_title(title)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save if path provided
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig

def plot_fitted_trajectories(
    x_values: np.ndarray,
    y_values: np.ndarray,
    fitted_values: np.ndarray,
    uncertainty: Optional[np.ndarray] = None,
    gene_name: Optional[str] = None,
    gene_idx: Optional[int] = None,
    model_name: str = "Fitted Model",
    figsize: Tuple[int, int] = (10, 6),
    colors: Tuple[str, str, str] = ('#3498db', '#e74c3c', '#3498db'),
    alpha_raw: float = 0.5,
    alpha_ci: float = 0.2,
    linewidth: float = 2.0,
    title: Optional[str] = None,
    legend_loc: str = 'best',
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    ax: Optional[plt.Axes] = None
) -> plt.Figure:
    """
    Plot fitted trajectories with original data points.
    
    Parameters
    ----------
    x_values : numpy.ndarray
        X values (time points)
    y_values : numpy.ndarray
        Y values (original data)
    fitted_values : numpy.ndarray
        Fitted Y values
    uncertainty : numpy.ndarray, optional
        Uncertainty bounds (std or confidence intervals)
    gene_name : str, optional
        Name of the gene
    gene_idx : int, optional
        Index of the gene
    model_name : str, optional
        Name of the fitted model
    figsize : tuple, optional
        Figure size (width, height)
    colors : tuple, optional
        Colors for (raw data, fitted curve, confidence interval)
    alpha_raw : float, optional
        Alpha for raw data points
    alpha_ci : float, optional
        Alpha for confidence interval
    linewidth : float, optional
        Line width for fitted curve
    title : str, optional
        Plot title
    legend_loc : str, optional
        Legend location
    xlim : tuple, optional
        X-axis limits
    ylim : tuple, optional
        Y-axis limits
    save_path : str, optional
        Path to save figure
    ax : matplotlib.axes.Axes, optional
        Axes to plot on (if None, creates new figure)
        
    Returns
    -------
    matplotlib.figure.Figure
        The figure with plotted trajectories
    """
    # Create figure if no axes provided
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    
    # Plot original data points
    ax.scatter(x_values, y_values, color=colors[0], alpha=alpha_raw, s=30, label='Original Data')
    
    # Plot fitted curve
    ax.plot(x_values, fitted_values, color=colors[1], linewidth=linewidth, label=model_name)
    
    # Plot uncertainty if provided
    if uncertainty is not None:
        if len(uncertainty.shape) == 1:
            # Single array of uncertainties (like standard deviation)
            upper = fitted_values + uncertainty
            lower = fitted_values - uncertainty
        elif len(uncertainty.shape) == 2 and uncertainty.shape[0] == 2:
            # Upper and lower bounds provided
            lower, upper = uncertainty
        else:
            raise ValueError(f"Unexpected uncertainty shape: {uncertainty.shape}. Expected (n,) or (2, n).")
        
        ax.fill_between(x_values, lower, upper, color=colors[2], alpha=alpha_ci, label='95% CI')
    
    # Set title
    if title is None:
        if gene_name:
            title = f"Fitted Trajectory for {gene_name}"
        elif gene_idx is not None:
            title = f"Fitted Trajectory for Gene {gene_idx}"
        else:
            title = "Fitted Trajectory"
    ax.set_title(title)
    
    # Set labels
    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("Expression")
    
    # Set limits if provided
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    
    # Add legend
    ax.legend(loc=legend_loc)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Tight layout
    plt.tight_layout()
    
    # Save if path provided
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig

def plot_fitting_comparison(
    x_values: np.ndarray,
    y_values: np.ndarray,
    fitted_models: Dict[str, np.ndarray],
    uncertainties: Optional[Dict[str, np.ndarray]] = None,
    metrics: Optional[Dict[str, Dict[str, float]]] = None,
    gene_name: Optional[str] = None,
    gene_idx: Optional[int] = None,
    figsize: Tuple[int, int] = (12, 8),
    colors: Optional[List[str]] = None,
    alpha_raw: float = 0.5,
    alpha_ci: float = 0.2,
    linewidth: float = 2.0,
    title: Optional[str] = None,
    add_metric_table: bool = True,
    table_metrics: Optional[List[str]] = None,
    legend_loc: str = 'best',
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Plot comparison of multiple fitted models.
    
    Parameters
    ----------
    x_values : numpy.ndarray
        X values (time points)
    y_values : numpy.ndarray
        Y values (original data)
    fitted_models : dict
        Dictionary of model_name: fitted_values
    uncertainties : dict, optional
        Dictionary of model_name: uncertainty bounds
    metrics : dict, optional
        Dictionary of model_name: {metric_name: value}
    gene_name : str, optional
        Name of the gene
    gene_idx : int, optional
        Index of the gene
    figsize : tuple, optional
        Figure size (width, height)
    colors : list, optional
        List of colors for models
    alpha_raw : float, optional
        Alpha for raw data points
    alpha_ci : float, optional
        Alpha for confidence intervals
    linewidth : float, optional
        Line width for fitted curves
    title : str, optional
        Plot title
    add_metric_table : bool, optional
        Whether to add metrics table
    table_metrics : list, optional
        List of metrics to include in table
    legend_loc : str, optional
        Legend location
    xlim : tuple, optional
        X-axis limits
    ylim : tuple, optional
        Y-axis limits
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    matplotlib.figure.Figure
        The figure with comparison plots
    """
    # Validate inputs
    if uncertainties is None:
        uncertainties = {}
    
    if metrics is None:
        metrics = {}
        
    if table_metrics is None and metrics:
        # Use all available metrics from the first model
        first_model = list(metrics.keys())[0]
        table_metrics = list(metrics[first_model].keys())
    
    # Set up colors
    if colors is None:
        # Use default color cycle
        color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
        colors = [color_cycle[i % len(color_cycle)] for i in range(len(fitted_models))]
    
    # Determine figure layout based on whether we need a metrics table
    if add_metric_table and metrics and table_metrics:
        # Figure with main plot and table
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(1, 2, width_ratios=[3, 1])
        ax_main = fig.add_subplot(gs[0, 0])
        ax_table = fig.add_subplot(gs[0, 1])
    else:
        # Just the main plot
        fig, ax_main = plt.subplots(figsize=figsize)
    
    # Plot original data points
    ax_main.scatter(x_values, y_values, color='gray', alpha=alpha_raw, s=30, label='Original Data')
    
    # Plot each model
    for i, (model_name, fitted_values) in enumerate(fitted_models.items()):
        color = colors[i % len(colors)]
        
        # Plot fitted curve
        ax_main.plot(x_values, fitted_values, color=color, linewidth=linewidth, label=model_name)
        
        # Plot uncertainty if available
        if model_name in uncertainties:
            uncertainty = uncertainties[model_name]
            
            if len(uncertainty.shape) == 1:
                # Single array of uncertainties (like standard deviation)
                upper = fitted_values + uncertainty
                lower = fitted_values - uncertainty
            elif len(uncertainty.shape) == 2 and uncertainty.shape[0] == 2:
                # Upper and lower bounds provided
                lower, upper = uncertainty
            else:
                continue  # Skip if format is not recognized
            
            ax_main.fill_between(x_values, lower, upper, color=color, alpha=alpha_ci)
    
    # Set title
    if title is None:
        if gene_name:
            title = f"Model Comparison for {gene_name}"
        elif gene_idx is not None:
            title = f"Model Comparison for Gene {gene_idx}"
        else:
            title = "Model Comparison"
    ax_main.set_title(title)
    
    # Set labels
    ax_main.set_xlabel("Pseudotime")
    ax_main.set_ylabel("Expression")
    
    # Set limits if provided
    if xlim is not None:
        ax_main.set_xlim(xlim)
    if ylim is not None:
        ax_main.set_ylim(ylim)
    
    # Add legend
    ax_main.legend(loc=legend_loc)
    
    # Add grid
    ax_main.grid(True, alpha=0.3)
    
    # Add metrics table if requested
    if add_metric_table and metrics and table_metrics and 'ax_table' in locals():
        # Create data for table
        table_data = []
        model_names = list(fitted_models.keys())
        
        # Add header row
        header = ['Metric'] + model_names
        table_data.append(header)
        
        # Add metric rows
        for metric in table_metrics:
            row = [metric]
            for model in model_names:
                if model in metrics and metric in metrics[model]:
                    value = metrics[model][metric]
                    # Format based on value type
                    if isinstance(value, float):
                        formatted = f"{value:.4f}"
                    else:
                        formatted = str(value)
                    row.append(formatted)
                else:
                    row.append('-')
            table_data.append(row)
        
        # Create table
        table = ax_table.table(
            cellText=table_data[1:],
            colLabels=table_data[0],
            loc='center',
            cellLoc='center'
        )
        
        # Style table
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        
        # Hide axis
        ax_table.axis('off')
        
        # Add title to table
        ax_table.set_title("Metrics Comparison")
    
    # Tight layout
    plt.tight_layout()
    
    # Save if path provided
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig

def visualize_fitting_results(
    standard_results: Dict[str, Any],
    optimized_results: Dict[str, Any],
    top_genes_data: List[np.ndarray],
    top_gene_names: List[str],
    time_points: np.ndarray,
    output_dir: Union[str, Path],
    max_genes_to_plot: int = 10,
    figsize: Tuple[int, int] = (12, 8)
) -> Dict[str, str]:
    """
    Visualize fitting results comparing standard and DTW-optimized approaches.
    
    Parameters
    ----------
    standard_results : dict
        Results from standard fitting approach
    optimized_results : dict
        Results from DTW-optimized fitting approach
    top_genes_data : list
        List of arrays containing gene data
    top_gene_names : list
        List of gene names
    time_points : numpy.ndarray
        Time points for the data
    output_dir : str or Path
        Directory to save visualizations
    max_genes_to_plot : int, optional
        Maximum number of genes to plot
    figsize : tuple, optional
        Figure size (width, height)
        
    Returns
    -------
    dict
        Dictionary with paths to visualization files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Limit number of genes to plot
    n_genes_to_plot = min(len(top_gene_names), max_genes_to_plot)
    
    # Create mapping for paths
    visualization_paths = {}
    
    # Plot comparison for each gene
    for i in range(n_genes_to_plot):
        gene_name = top_gene_names[i]
        gene_data = top_genes_data[i]
        
        # Calculate mean trajectory for this gene
        mean_trajectory = np.mean(gene_data, axis=0)[:, 0]
        
        # Extract fitted trajectories for this gene
        standard_traj = standard_results['fitted_trajectories'][:, i]
        optimized_traj = optimized_results['fitted_trajectories'][:, i]
        
        # Get metrics
        standard_dtw = standard_results['dtw_distances'][i]
        optimized_dtw = optimized_results['dtw_distances'][i]
        standard_smoothing = standard_results['smoothing_values'][i]
        optimized_smoothing = optimized_results['smoothing_values'][i]
        
        # Create metrics dictionary
        metrics = {
            'Standard': {
                'DTW Distance': standard_dtw,
                'Smoothing': standard_smoothing
            },
            'DTW-Optimized': {
                'DTW Distance': optimized_dtw,
                'Smoothing': optimized_smoothing,
                'Improvement (%)': 100 * (standard_dtw - optimized_dtw) / standard_dtw if standard_dtw > 0 else 0
            }
        }
        
        # Create model dictionary
        fitted_models = {
            'Standard': standard_traj,
            'DTW-Optimized': optimized_traj
        }
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        gs = plt.GridSpec(2, 1, height_ratios=[3, 1])
        ax_main = fig.add_subplot(gs[0])
        ax_stats = fig.add_subplot(gs[1])
        
        # Plot data points
        for j in range(min(gene_data.shape[0], 10)):  # Plot max 10 samples for clarity
            ax_main.plot(time_points, gene_data[j, :, 0], 'o', alpha=0.3, markersize=3)
        
        # Plot mean trajectory
        ax_main.plot(time_points, mean_trajectory, 'k--', alpha=0.7, label='Mean')
        
        # Plot fitted models
        ax_main.plot(standard_results['time_points'], standard_traj, 'b-', linewidth=2, label='Standard')
        ax_main.plot(optimized_results['time_points'], optimized_traj, 'r-', linewidth=2, label='DTW-Optimized')
        
        # Set title and labels
        ax_main.set_title(f"Trajectory Fitting Comparison for {gene_name}")
        ax_main.set_xlabel("Pseudotime")
        ax_main.set_ylabel("Expression")
        ax_main.legend()
        ax_main.grid(alpha=0.3)
        
        # Create table for metrics
        table_data = [
            ['Method', 'DTW Distance', 'Smoothing', 'Improvement (%)'],
            ['Standard', f"{standard_dtw:.4f}", f"{standard_smoothing:.4f}", "-"],
            ['DTW-Optimized', f"{optimized_dtw:.4f}", f"{optimized_smoothing:.4f}", 
             f"{100 * (standard_dtw - optimized_dtw) / standard_dtw:.2f}%" if standard_dtw > 0 else "N/A"]
        ]
        
        # Hide axes for table subplot
        ax_stats.axis('off')
        
        # Create table
        table = ax_stats.table(
            cellText=table_data[1:],
            colLabels=table_data[0],
            loc='center',
            cellLoc='center'
        )
        
        # Style table
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        
        # Save figure
        save_path = output_dir / f"fitting_comparison_{gene_name}.png"
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        # Store path
        visualization_paths[gene_name] = str(save_path)
    
    # Create summary figure with DTW distances for all genes
    fig, ax = plt.subplots(figsize=figsize)
    
    n_genes = min(len(top_gene_names), 20)  # Limit to 20 genes for readability
    indices = np.arange(n_genes)
    
    # Plot bar chart
    width = 0.35
    ax.bar(indices - width/2, standard_results['dtw_distances'][:n_genes], width, label='Standard')
    ax.bar(indices + width/2, optimized_results['dtw_distances'][:n_genes], width, label='DTW-Optimized')
    
    # Add labels and title
    ax.set_xlabel('Gene')
    ax.set_ylabel('DTW Distance')
    ax.set_title('Comparison of DTW Distances by Gene')
    ax.set_xticks(indices)
    ax.set_xticklabels(top_gene_names[:n_genes], rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Save figure
    summary_path = output_dir / "dtw_distance_comparison.png"
    plt.tight_layout()
    plt.savefig(summary_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Add summary to paths
    visualization_paths['summary'] = str(summary_path)
    
    # Create smoothing values comparison
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot bar chart
    ax.bar(indices - width/2, standard_results['smoothing_values'][:n_genes], width, label='Standard')
    ax.bar(indices + width/2, optimized_results['smoothing_values'][:n_genes], width, label='DTW-Optimized')
    
    # Add labels and title
    ax.set_xlabel('Gene')
    ax.set_ylabel('Smoothing Value')
    ax.set_title('Comparison of Smoothing Values by Gene')
    ax.set_xticks(indices)
    ax.set_xticklabels(top_gene_names[:n_genes], rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Save figure
    smoothing_path = output_dir / "smoothing_comparison.png"
    plt.tight_layout()
    plt.savefig(smoothing_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Add smoothing to paths
    visualization_paths['smoothing'] = str(smoothing_path)
    
    return visualization_paths

def create_fitting_summary(
    standard_results: Dict[str, Any],
    optimized_results: Dict[str, Any],
    top_gene_names: List[str],
    top_genes_data: List[np.ndarray],
    output_file: Union[str, Path],
    adata_shape: Optional[Tuple[int, int]] = None,
    reshaped_data_shape: Optional[Tuple[int, int, int]] = None,
    batch_names: Optional[List[str]] = None
) -> str:
    """
    Create a summary report of fitting results.
    
    Parameters
    ----------
    standard_results : dict
        Results from standard fitting approach
    optimized_results : dict
        Results from DTW-optimized fitting approach
    top_gene_names : list
        List of gene names
    top_genes_data : list
        List of arrays containing gene data
    output_file : str or Path
        Path to save summary report
    adata_shape : tuple, optional
        Shape of original AnnData object (n_cells, n_genes)
    reshaped_data_shape : tuple, optional
        Shape of reshaped 3D data (n_batches, n_timepoints, n_genes)
    batch_names : list, optional
        Names of batches
        
    Returns
    -------
    str
        Path to summary report
    """
    output_file = Path(output_file)
    output_dir = output_file.parent
    output_dir.mkdir(exist_ok=True, parents=True)
    
    try:
        with open(output_file, 'w') as f:
            # Write header
            f.write("=" * 80 + "\n")
            f.write("TRAJECTORY FITTING ANALYSIS SUMMARY\n")
            f.write("=" * 80 + "\n\n")
            
            # Write date and time
            from datetime import datetime
            f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Write dataset information
            f.write("DATASET INFORMATION\n")
            f.write("-" * 80 + "\n")
            
            if adata_shape is not None:
                f.write(f"Original data shape: {adata_shape[0]} cells × {adata_shape[1]} genes\n")
            
            if reshaped_data_shape is not None:
                f.write(f"Reshaped data: {reshaped_data_shape[0]} batches × {reshaped_data_shape[1]} timepoints × {reshaped_data_shape[2]} genes\n")
            
            if batch_names is not None:
                f.write(f"Batches: {', '.join(batch_names)}\n")
            
            f.write(f"Top genes analyzed: {len(top_gene_names)}\n")
            f.write("\n")
            
            # Write overall results
            f.write("OVERALL RESULTS\n")
            f.write("-" * 80 + "\n")
            
            # Calculate mean DTW distances
            standard_mean_dtw = np.mean(standard_results['dtw_distances'])
            optimized_mean_dtw = np.mean(optimized_results['dtw_distances'])
            
            f.write(f"Standard approach - mean DTW distance: {standard_mean_dtw:.4f}\n")
            f.write(f"DTW-optimized approach - mean DTW distance: {optimized_mean_dtw:.4f}\n")
            
            improvement = standard_mean_dtw - optimized_mean_dtw
            percent_improvement = 100 * improvement / standard_mean_dtw if standard_mean_dtw > 0 else 0
            
            f.write(f"Absolute improvement: {improvement:.4f}\n")
            f.write(f"Percentage improvement: {percent_improvement:.2f}%\n\n")
            
            # Write detailed gene results
            f.write("GENE-LEVEL RESULTS\n")
            f.write("-" * 80 + "\n")
            f.write("Gene\tStandard DTW\tOptimized DTW\tImprovement\tStd Smoothing\tOpt Smoothing\n")
            
            for i, gene_name in enumerate(top_gene_names):
                standard_dtw = standard_results['dtw_distances'][i]
                optimized_dtw = optimized_results['dtw_distances'][i]
                improvement = standard_dtw - optimized_dtw
                percent_imp = 100 * improvement / standard_dtw if standard_dtw > 0 else 0
                
                standard_smoothing = standard_results['smoothing_values'][i]
                optimized_smoothing = optimized_results['smoothing_values'][i]
                
                f.write(f"{gene_name}\t{standard_dtw:.4f}\t{optimized_dtw:.4f}\t{percent_imp:.2f}%\t{standard_smoothing:.4f}\t{optimized_smoothing:.4f}\n")
            
            f.write("\n")
            
            # Write genes with largest improvement
            improvements = [(i, standard_results['dtw_distances'][i] - optimized_results['dtw_distances'][i]) 
                           for i in range(len(top_gene_names))]
            improvements.sort(key=lambda x: x[1], reverse=True)
            
            f.write("TOP 5 GENES WITH LARGEST IMPROVEMENT\n")
            f.write("-" * 80 + "\n")
            f.write("Gene\tStandard DTW\tOptimized DTW\tImprovement\tPercent Improvement\n")
            
            for i, imp in improvements[:5]:
                gene_name = top_gene_names[i]
                standard_dtw = standard_results['dtw_distances'][i]
                optimized_dtw = optimized_results['dtw_distances'][i]
                percent_imp = 100 * imp / standard_dtw if standard_dtw > 0 else 0
                
                f.write(f"{gene_name}\t{standard_dtw:.4f}\t{optimized_dtw:.4f}\t{imp:.4f}\t{percent_imp:.2f}%\n")
                
            # Write conclusions
            f.write("\n")
            f.write("CONCLUSIONS\n")
            f.write("-" * 80 + "\n")
            
            if percent_improvement > 10:
                f.write("The DTW optimization significantly improved the fitting quality.\n")
            elif percent_improvement > 0:
                f.write("The DTW optimization slightly improved the fitting quality.\n")
            else:
                f.write("The DTW optimization did not improve the overall fitting quality.\n")
                
            # Add recommendation based on results
            avg_standard_smoothing = np.mean(standard_results['smoothing_values'])
            avg_optimized_smoothing = np.mean(optimized_results['smoothing_values'])
            
            if avg_optimized_smoothing < avg_standard_smoothing:
                f.write("Recommended smoother fits (lower smoothing values) for this dataset.\n")
            else:
                f.write("Recommended less smooth fits (higher smoothing values) for this dataset.\n")
                
            # Add file information
            f.write("\n")
            f.write("FILE INFORMATION\n")
            f.write("-" * 80 + "\n")
            f.write(f"Summary report: {output_file}\n")
            
        return str(output_file)
        
    except Exception as e:
        # If there's an error, write to an error file
        error_file = output_dir / "fitting_summary_error.txt"
        with open(error_file, 'w') as f:
            f.write(f"Error creating summary: {str(e)}\n")
        return str(error_file) 