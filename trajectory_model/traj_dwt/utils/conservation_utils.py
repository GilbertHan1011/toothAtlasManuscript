import numpy as np
import warnings
from fastdtw import fastdtw
from scipy.spatial.distance import euclidean

def custom_distance(x, y):
    """
    Custom distance function for DTW that handles scalar values.
    
    Parameters
    ----------
    x : any
        First point
    y : any
        Second point
        
    Returns
    -------
    float
        Absolute difference between points
    """
    return float(abs(x - y))

def compute_dtw(x: np.ndarray, y: np.ndarray, radius: int = 3) -> float:
    """
    Compute DTW distance using the best available method.
    
    Parameters
    ----------
    x : np.ndarray
        First sequence
    y : np.ndarray
        Second sequence
    radius : int
        Radius parameter for fastdtw
        
    Returns
    -------
    float
        DTW distance
    """
    # Ensure input arrays are 1D
    x = np.ravel(x)
    y = np.ravel(y)
    
    # Use fastdtw if available (fastest)
    if FASTDTW_AVAILABLE:
        try:
            # Try with custom distance function
            distance, _ = fastdtw(x, y, dist=custom_distance, radius=radius)
            return distance
        except Exception as e:
            warnings.warn(f"Error using fastdtw: {e}. Falling back to custom implementation.")
            return _optimized_dtw(x, y, window=radius)
    
    # Try to use scipy's DTW if available
    try:
        from scipy.spatial.distance import dtw
        distance, _, _, _ = dtw(x, y, dist=euclidean)
        return distance
    except (ImportError, AttributeError):
        # If scipy's DTW is not available, use our optimized implementation
        return _optimized_dtw(x, y, window=radius)