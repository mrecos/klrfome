# KLRfome Python/JAX Refactor Specification

## Project Overview

Refactor the KLRfome (Kernel Logistic Regression on Focal Mean Embeddings) R package into a modern, high-performance Python library using JAX for GPU acceleration and automatic differentiation.

**Original R Package**: https://github.com/mrecos/klrfome
**Original Paper**: Harris, M.D. (2019). KLRfome - Kernel Logistic Regression on Focal Mean Embeddings. Journal of Open Source Software, 4(35), 722. https://doi.org/10.21105/joss.00722

---

## What KLRfome Does (Conceptual Background)

KLRfome solves a **Distribution Regression** problem for geographic/spatial prediction. The core insight is:

1. **The Problem**: We want to predict a binary outcome (e.g., site present/absent) across a landscape. However, each "site" is not a single point—it's an *area* characterized by a *distribution* of environmental measurements (elevation, slope, aspect, soil type, etc.).

2. **Traditional Approaches Fail Because**:
   - Collapsing each site to its centroid/mean loses distributional information
   - Treating each raster cell as independent violates spatial autocorrelation assumptions and inflates sample sizes

3. **The KLRfome Solution**:
   - Represent each site as a distribution of feature vectors (samples from within the site boundary)
   - Compute similarity between sites using kernel methods on these distributions (mean embeddings in RKHS)
   - Fit a Kernel Logistic Regression model on the similarity matrix
   - Predict using a focal (moving window) approach that computes similarity between each landscape neighborhood and the training data

---

## Core Algorithm Components

### 1. Data Structure

**Input Data Format**:
- Environmental raster stack (multiple co-registered rasters representing different variables)
- Site locations (polygons or point buffers defining "presence" locations)
- Background samples (random or stratified samples from non-site areas)

**Internal Data Structure**:
For each site/background location, we store multiple samples:
```
{
  "site_1": [[var1, var2, ..., varN], [var1, var2, ..., varN], ...],  # M samples
  "site_2": [[var1, var2, ..., varN], [var1, var2, ..., varN], ...],
  ...
  "background_1": [[var1, var2, ..., varN], ...],
  ...
}
```

Each site/background is represented by a *set* of feature vectors, not a single vector.

### 2. Similarity Kernel (Mean Embeddings)

The similarity between two distributions P and Q is computed via mean embeddings in a Reproducing Kernel Hilbert Space (RKHS):

```
K_distribution(P, Q) = k(μ_P, μ_Q)
```

Where:
- μ_P = (1/|P|) Σ_{x∈P} φ(x) is the mean embedding of distribution P
- φ(x) is the feature map induced by a base kernel k
- For Gaussian/RBF kernel: k(x, y) = exp(-||x - y||² / (2σ²))

**In practice**, for RBF kernels, the similarity between two sets of samples X = {x_1, ..., x_m} and Y = {y_1, ..., y_n} is:

```
K(X, Y) = (1/mn) Σ_i Σ_j k(x_i, y_j)
```

This is computationally expensive (O(mn) kernel evaluations per pair), which is why we need approximations.

### 3. Kernel Logistic Regression (KLR)

Given the N×N similarity matrix K between all training sites/backgrounds, fit a logistic regression in kernel space using Iteratively Reweighted Least Squares (IRLS):

**Model**:
```
P(y=1 | x) = 1 / (1 + exp(-f(x)))
f(x) = Σ_j α_j K(x, x_j)
```

**IRLS Algorithm**:
```
Initialize: α = 0, converged = False
While not converged:
    η = K @ α                           # Linear predictor
    p = 1 / (1 + exp(-η))               # Probabilities
    W = diag(p * (1 - p))               # Weights
    z = η + (y - p) / (p * (1 - p))     # Working response
    α_new = solve((K @ W @ K + λI), K @ W @ z)  # Ridge-regularized update
    Check convergence: ||α_new - α|| < tol
    α = α_new
```

**Hyperparameters**:
- `sigma`: Bandwidth of the Gaussian/RBF kernel
- `lambda`: L2 regularization strength

### 4. Focal Window Prediction

Prediction is NOT done cell-by-cell. Instead, for each location in the prediction raster:

1. Extract an N×N window of cells around the target location
2. Treat this window as a "distribution" (set of feature vectors)
3. Compute similarity between this window distribution and all training sites/backgrounds
4. Apply the learned KLR model to get P(site | window)

**Hyperparameters**:
- `ngb`: Neighborhood size (e.g., 3 means 3×3 window, 5 means 5×5, etc.)

This focal approach:
- Captures the distributional nature of landforms
- Produces smoother, more realistic prediction surfaces
- Is embarrassingly parallel (each window is independent)

---

## Python/JAX Implementation Requirements

### Package Structure

```
klrfome/
├── __init__.py
├── data/
│   ├── __init__.py
│   ├── formats.py          # Data structure definitions and conversions
│   ├── sampling.py         # Extract samples from rasters at locations
│   └── simulation.py       # Generate simulated test data
├── kernels/
│   ├── __init__.py
│   ├── base.py             # Base kernel classes/protocols
│   ├── rbf.py              # Exact RBF kernel
│   ├── rff.py              # Random Fourier Features approximation
│   └── distribution.py     # Distribution-level kernels (mean embedding)
├── models/
│   ├── __init__.py
│   ├── klr.py              # Kernel Logistic Regression
│   └── utils.py            # IRLS solver, convergence checking
├── prediction/
│   ├── __init__.py
│   ├── focal.py            # Focal window prediction
│   └── parallel.py         # Parallelization utilities
├── visualization/
│   ├── __init__.py
│   └── plots.py            # Plotting utilities
├── io/
│   ├── __init__.py
│   ├── raster.py           # Raster I/O (rasterio integration)
│   └── vector.py           # Vector I/O (geopandas integration)
└── utils/
    ├── __init__.py
    ├── validation.py       # Cross-validation, metrics
    └── scaling.py          # Feature scaling/normalization
```

### Core Dependencies

```
# Core computation
jax>=0.4.0
jaxlib>=0.4.0
numpy>=1.24.0
scipy>=1.10.0

# Geospatial
rasterio>=1.3.0
geopandas>=0.12.0
shapely>=2.0.0
pyproj>=3.4.0

# Data handling
pandas>=2.0.0
xarray>=2023.0.0  # Optional but recommended for raster stacks

# Visualization
matplotlib>=3.7.0
seaborn>=0.12.0

# Utilities
tqdm>=4.65.0
typing_extensions>=4.5.0

# Testing
pytest>=7.0.0
pytest-cov>=4.0.0
hypothesis>=6.0.0  # Property-based testing
```

---

## Detailed Implementation Specifications

### 1. Data Structures (`klrfome/data/formats.py`)

```python
from dataclasses import dataclass
from typing import Dict, List, Optional, Literal
import jax.numpy as jnp
from jaxtyping import Array, Float

@dataclass
class SampleCollection:
    """
    Represents samples from a single site or background location.
    
    Attributes:
        samples: Array of shape (n_samples, n_features)
        label: 1 for site, 0 for background
        id: Unique identifier for this location
        metadata: Optional dict for additional info (coordinates, etc.)
    """
    samples: Float[Array, "n_samples n_features"]
    label: Literal[0, 1]
    id: str
    metadata: Optional[Dict] = None
    
    @property
    def n_samples(self) -> int:
        return self.samples.shape[0]
    
    @property
    def n_features(self) -> int:
        return self.samples.shape[1]
    
    def mean_embedding(self) -> Float[Array, "n_features"]:
        """Compute the empirical mean of samples."""
        return jnp.mean(self.samples, axis=0)


@dataclass
class TrainingData:
    """
    Complete training dataset for KLRfome.
    
    Attributes:
        collections: List of SampleCollection objects
        feature_names: Names of the environmental variables
        crs: Coordinate reference system (from rasterio/pyproj)
    """
    collections: List[SampleCollection]
    feature_names: List[str]
    crs: Optional[str] = None
    
    @property
    def n_locations(self) -> int:
        return len(self.collections)
    
    @property
    def labels(self) -> Float[Array, "n_locations"]:
        return jnp.array([c.label for c in self.collections])
    
    @property
    def n_sites(self) -> int:
        return sum(1 for c in self.collections if c.label == 1)
    
    @property
    def n_background(self) -> int:
        return sum(1 for c in self.collections if c.label == 0)
    
    def train_test_split(
        self, 
        test_fraction: float = 0.2, 
        stratify: bool = True,
        seed: int = 42
    ) -> tuple['TrainingData', 'TrainingData']:
        """Split into training and testing sets."""
        ...


@dataclass 
class RasterStack:
    """
    Wrapper around a stack of co-registered rasters.
    
    Provides convenient access for focal window extraction.
    """
    data: Float[Array, "n_bands height width"]
    transform: 'rasterio.Affine'
    crs: str
    band_names: List[str]
    nodata: Optional[float] = None
    
    @classmethod
    def from_files(cls, file_paths: List[str]) -> 'RasterStack':
        """Load from multiple single-band rasters."""
        ...
    
    @classmethod
    def from_multiband(cls, file_path: str) -> 'RasterStack':
        """Load from a single multi-band raster."""
        ...
    
    def extract_window(
        self, 
        row: int, 
        col: int, 
        size: int
    ) -> Float[Array, "size size n_bands"]:
        """Extract an NxN window centered at (row, col)."""
        ...
    
    def extract_at_points(
        self,
        points: 'geopandas.GeoDataFrame',
        buffer_radius: Optional[float] = None,
        n_samples: int = 10
    ) -> List[SampleCollection]:
        """Extract samples at point locations, optionally with buffer."""
        ...
```

### 2. Kernel Implementations (`klrfome/kernels/`)

#### Base Kernel Protocol (`base.py`)

```python
from typing import Protocol, runtime_checkable
from jaxtyping import Array, Float

@runtime_checkable
class Kernel(Protocol):
    """Protocol for kernel functions."""
    
    def __call__(
        self, 
        X: Float[Array, "n d"], 
        Y: Float[Array, "m d"]
    ) -> Float[Array, "n m"]:
        """Compute kernel matrix between X and Y."""
        ...
    
    @property
    def sigma(self) -> float:
        """Kernel bandwidth parameter."""
        ...


@runtime_checkable
class ApproximateKernel(Protocol):
    """Protocol for kernels with explicit feature maps."""
    
    def feature_map(
        self, 
        X: Float[Array, "n d"]
    ) -> Float[Array, "n D"]:
        """Map inputs to approximate feature space."""
        ...
```

#### Exact RBF Kernel (`rbf.py`)

```python
import jax.numpy as jnp
from jax import jit
from functools import partial

class RBFKernel:
    """
    Radial Basis Function (Gaussian) kernel.
    
    k(x, y) = exp(-||x - y||² / (2σ²))
    
    Parameters:
        sigma: Bandwidth parameter (length scale)
    """
    
    def __init__(self, sigma: float = 1.0):
        self._sigma = sigma
    
    @property
    def sigma(self) -> float:
        return self._sigma
    
    @partial(jit, static_argnums=(0,))
    def __call__(
        self, 
        X: Float[Array, "n d"], 
        Y: Float[Array, "m d"]
    ) -> Float[Array, "n m"]:
        """
        Compute RBF kernel matrix.
        
        Uses the identity:
        ||x - y||² = ||x||² + ||y||² - 2<x, y>
        """
        X_sqnorm = jnp.sum(X ** 2, axis=1, keepdims=True)
        Y_sqnorm = jnp.sum(Y ** 2, axis=1, keepdims=True)
        sq_distances = X_sqnorm + Y_sqnorm.T - 2 * jnp.dot(X, Y.T)
        return jnp.exp(-sq_distances / (2 * self._sigma ** 2))
    
    def diagonal(self, X: Float[Array, "n d"]) -> Float[Array, "n"]:
        """Diagonal of K(X, X) - always 1 for RBF."""
        return jnp.ones(X.shape[0])
```

#### Random Fourier Features (`rff.py`)

```python
import jax.numpy as jnp
import jax.random as random
from jax import jit
from functools import partial

class RandomFourierFeatures:
    """
    Random Fourier Features approximation to RBF kernel.
    
    Approximates: k(x, y) ≈ φ(x)ᵀφ(y)
    
    Where φ(x) = sqrt(2/D) * cos(Wx + b)
    
    This transforms the kernel computation from O(n²) to O(nD) where
    D is the number of random features.
    
    Reference:
        Rahimi & Recht (2007). "Random Features for Large-Scale Kernel Machines"
    
    Parameters:
        sigma: Bandwidth of the RBF kernel to approximate
        n_features: Number of random features (D). Higher = better approximation.
        seed: Random seed for reproducibility
    """
    
    def __init__(
        self, 
        sigma: float = 1.0, 
        n_features: int = 256,
        seed: int = 42
    ):
        self._sigma = sigma
        self._n_features = n_features
        self._seed = seed
        self._W = None  # Lazily initialized
        self._b = None
        self._input_dim = None
    
    @property
    def sigma(self) -> float:
        return self._sigma
    
    @property
    def n_features(self) -> int:
        return self._n_features
    
    def _initialize_weights(self, input_dim: int):
        """Initialize random weights for the feature map."""
        if self._W is not None and self._input_dim == input_dim:
            return
        
        key = random.PRNGKey(self._seed)
        key_W, key_b = random.split(key)
        
        # W ~ N(0, 1/σ²) for RBF kernel
        self._W = random.normal(key_W, (input_dim, self._n_features)) / self._sigma
        # b ~ Uniform(0, 2π)
        self._b = random.uniform(key_b, (self._n_features,), minval=0, maxval=2*jnp.pi)
        self._input_dim = input_dim
    
    def feature_map(
        self, 
        X: Float[Array, "n d"]
    ) -> Float[Array, "n D"]:
        """
        Compute random Fourier features.
        
        φ(x) = sqrt(2/D) * cos(Wx + b)
        """
        self._initialize_weights(X.shape[1])
        projection = jnp.dot(X, self._W) + self._b
        return jnp.sqrt(2.0 / self._n_features) * jnp.cos(projection)
    
    def __call__(
        self, 
        X: Float[Array, "n d"], 
        Y: Float[Array, "m d"]
    ) -> Float[Array, "n m"]:
        """
        Approximate kernel matrix via random features.
        
        K(X, Y) ≈ φ(X) @ φ(Y).T
        """
        phi_X = self.feature_map(X)
        phi_Y = self.feature_map(Y)
        return jnp.dot(phi_X, phi_Y.T)
    
    def self_similarity(
        self, 
        X: Float[Array, "n d"]
    ) -> Float[Array, "n n"]:
        """Compute K(X, X) efficiently."""
        phi_X = self.feature_map(X)
        return jnp.dot(phi_X, phi_X.T)
```

#### Distribution Kernel (`distribution.py`)

```python
import jax.numpy as jnp
from jax import jit, vmap
from typing import Union

class MeanEmbeddingKernel:
    """
    Kernel on distributions via mean embeddings.
    
    Given two sets of samples X = {x_1, ..., x_m} and Y = {y_1, ..., y_n},
    computes similarity as:
    
    K(X, Y) = (1/mn) Σ_i Σ_j k(x_i, y_j)
    
    This is equivalent to the inner product of mean embeddings in RKHS:
    <μ_X, μ_Y>_H where μ_X = (1/m) Σ_i φ(x_i)
    
    Parameters:
        base_kernel: The point-level kernel (e.g., RBFKernel or RandomFourierFeatures)
        use_rff: If True and base_kernel supports it, use RFF for efficiency
    """
    
    def __init__(
        self, 
        base_kernel: Union[RBFKernel, RandomFourierFeatures],
    ):
        self.base_kernel = base_kernel
        self._use_rff = hasattr(base_kernel, 'feature_map')
    
    def __call__(
        self,
        X: Float[Array, "m d"],
        Y: Float[Array, "n d"]
    ) -> float:
        """
        Compute distribution similarity between sample sets X and Y.
        
        Returns a scalar similarity value.
        """
        if self._use_rff:
            # Efficient: compute mean embeddings in feature space
            phi_X = self.base_kernel.feature_map(X)  # (m, D)
            phi_Y = self.base_kernel.feature_map(Y)  # (n, D)
            mean_X = jnp.mean(phi_X, axis=0)  # (D,)
            mean_Y = jnp.mean(phi_Y, axis=0)  # (D,)
            return jnp.dot(mean_X, mean_Y)
        else:
            # Exact: compute full kernel matrix and average
            K = self.base_kernel(X, Y)  # (m, n)
            return jnp.mean(K)
    
    def build_similarity_matrix(
        self,
        collections: List[SampleCollection]
    ) -> Float[Array, "N N"]:
        """
        Build the N×N similarity matrix between all collections.
        
        This is the core computation for KLR fitting.
        """
        n = len(collections)
        
        if self._use_rff:
            # Efficient batch computation with RFF
            # First, compute mean embeddings for all collections
            mean_embeddings = []
            for coll in collections:
                phi = self.base_kernel.feature_map(coll.samples)
                mean_embeddings.append(jnp.mean(phi, axis=0))
            
            mean_embeddings = jnp.stack(mean_embeddings)  # (N, D)
            # Similarity matrix is just dot products of mean embeddings
            return jnp.dot(mean_embeddings, mean_embeddings.T)
        
        else:
            # Exact computation (slower)
            K = jnp.zeros((n, n))
            for i in range(n):
                for j in range(i, n):
                    k_ij = self(collections[i].samples, collections[j].samples)
                    K = K.at[i, j].set(k_ij)
                    K = K.at[j, i].set(k_ij)
            return K
```

### 3. Kernel Logistic Regression (`klrfome/models/klr.py`)

```python
import jax.numpy as jnp
from jax import jit, grad
from jax.scipy.linalg import solve
from typing import Tuple, Optional, NamedTuple
from dataclasses import dataclass
import warnings

class KLRFitResult(NamedTuple):
    """Result of KLR fitting."""
    alpha: Float[Array, "n"]  # Dual coefficients
    converged: bool
    n_iterations: int
    final_loss: float


@dataclass
class KernelLogisticRegression:
    """
    Kernel Logistic Regression with IRLS solver.
    
    Fits the model:
        P(y=1 | x) = σ(Σ_j α_j K(x, x_j))
    
    where σ is the sigmoid function and K is the kernel.
    
    Parameters:
        lambda_reg: L2 regularization strength
        max_iter: Maximum IRLS iterations
        tol: Convergence tolerance for alpha
        min_prob: Minimum probability to avoid numerical issues (clipping)
    """
    lambda_reg: float = 1.0
    max_iter: int = 100
    tol: float = 1e-6
    min_prob: float = 1e-7
    
    def fit(
        self,
        K: Float[Array, "n n"],
        y: Float[Array, "n"],
        alpha_init: Optional[Float[Array, "n"]] = None
    ) -> KLRFitResult:
        """
        Fit KLR model using IRLS.
        
        Parameters:
            K: Precomputed kernel/similarity matrix
            y: Binary labels (0 or 1)
            alpha_init: Initial alpha values (default: zeros)
        
        Returns:
            KLRFitResult with fitted coefficients and diagnostics
        """
        n = K.shape[0]
        alpha = alpha_init if alpha_init is not None else jnp.zeros(n)
        
        for iteration in range(self.max_iter):
            # Compute probabilities
            eta = K @ alpha
            prob = self._sigmoid(eta)
            prob = jnp.clip(prob, self.min_prob, 1 - self.min_prob)
            
            # IRLS weights
            W = prob * (1 - prob)
            
            # Working response
            z = eta + (y - prob) / W
            
            # Weighted least squares update with regularization
            # Solve: (K W K + λI) α = K W z
            KW = K * W[None, :]  # Broadcasting for diagonal W
            lhs = KW @ K + self.lambda_reg * jnp.eye(n)
            rhs = KW @ z
            
            alpha_new = solve(lhs, rhs, assume_a='pos')
            
            # Check convergence
            delta = jnp.max(jnp.abs(alpha_new - alpha))
            alpha = alpha_new
            
            if delta < self.tol:
                loss = self._compute_loss(K, y, alpha)
                return KLRFitResult(alpha, True, iteration + 1, loss)
        
        warnings.warn(f"KLR did not converge in {self.max_iter} iterations")
        loss = self._compute_loss(K, y, alpha)
        return KLRFitResult(alpha, False, self.max_iter, loss)
    
    def predict_proba(
        self,
        K_new: Float[Array, "m n"],
        alpha: Float[Array, "n"]
    ) -> Float[Array, "m"]:
        """
        Predict probabilities for new data.
        
        Parameters:
            K_new: Kernel matrix between new points and training points
                   Shape: (n_new, n_train)
            alpha: Fitted dual coefficients
        
        Returns:
            Predicted probabilities of class 1
        """
        eta = K_new @ alpha
        return self._sigmoid(eta)
    
    def predict(
        self,
        K_new: Float[Array, "m n"],
        alpha: Float[Array, "n"],
        threshold: float = 0.5
    ) -> Float[Array, "m"]:
        """Predict binary labels."""
        proba = self.predict_proba(K_new, alpha)
        return (proba >= threshold).astype(jnp.int32)
    
    @staticmethod
    def _sigmoid(x: Float[Array, "..."]) -> Float[Array, "..."]:
        """Numerically stable sigmoid."""
        return jnp.where(
            x >= 0,
            1 / (1 + jnp.exp(-x)),
            jnp.exp(x) / (1 + jnp.exp(x))
        )
    
    def _compute_loss(
        self,
        K: Float[Array, "n n"],
        y: Float[Array, "n"],
        alpha: Float[Array, "n"]
    ) -> float:
        """Compute regularized negative log-likelihood."""
        prob = self.predict_proba(K, alpha)
        prob = jnp.clip(prob, self.min_prob, 1 - self.min_prob)
        
        # Negative log-likelihood
        nll = -jnp.mean(y * jnp.log(prob) + (1 - y) * jnp.log(1 - prob))
        # Regularization
        reg = 0.5 * self.lambda_reg * alpha @ K @ alpha
        
        return nll + reg
```

### 4. Focal Window Prediction (`klrfome/prediction/focal.py`)

```python
import jax.numpy as jnp
from jax import jit, vmap, pmap
import jax.lax as lax
from functools import partial
from typing import Tuple, Optional
from tqdm import tqdm

@dataclass
class FocalPredictor:
    """
    Focal window prediction for KLRfome.
    
    Slides a window across the raster stack, computing similarity
    between each window and training data, then applying the fitted
    KLR model.
    
    Parameters:
        distribution_kernel: MeanEmbeddingKernel instance
        klr_model: Fitted KernelLogisticRegression
        alpha: Fitted coefficients from KLR
        training_collections: Original training SampleCollections
        window_size: Size of focal window (e.g., 3 for 3x3)
        use_gpu: Whether to use GPU acceleration
    """
    distribution_kernel: MeanEmbeddingKernel
    klr_alpha: Float[Array, "n_train"]
    training_data: TrainingData
    window_size: int = 3
    use_gpu: bool = True
    
    def __post_init__(self):
        # Precompute training mean embeddings for efficiency
        if self.distribution_kernel._use_rff:
            self._training_embeddings = self._compute_training_embeddings()
        else:
            self._training_embeddings = None
    
    def _compute_training_embeddings(self) -> Float[Array, "n_train D"]:
        """Precompute mean embeddings for all training collections."""
        embeddings = []
        for coll in self.training_data.collections:
            phi = self.distribution_kernel.base_kernel.feature_map(coll.samples)
            embeddings.append(jnp.mean(phi, axis=0))
        return jnp.stack(embeddings)
    
    def predict_window(
        self,
        window_samples: Float[Array, "m d"]
    ) -> float:
        """
        Predict probability for a single focal window.
        
        Parameters:
            window_samples: Samples from the focal window, shape (m, n_features)
                           where m = window_size * window_size (excluding nodata)
        
        Returns:
            Predicted probability of "site" class
        """
        if self._training_embeddings is not None:
            # RFF path: compute mean embedding and dot with training embeddings
            phi = self.distribution_kernel.base_kernel.feature_map(window_samples)
            window_embedding = jnp.mean(phi, axis=0)
            K_new = jnp.dot(window_embedding, self._training_embeddings.T)
        else:
            # Exact path: compute kernel with each training collection
            K_new = jnp.array([
                self.distribution_kernel(window_samples, coll.samples)
                for coll in self.training_data.collections
            ])
        
        # Apply KLR prediction
        eta = jnp.dot(K_new, self.klr_alpha)
        return 1 / (1 + jnp.exp(-eta))
    
    def predict_raster(
        self,
        raster_stack: RasterStack,
        batch_size: int = 1000,
        show_progress: bool = True
    ) -> Float[Array, "height width"]:
        """
        Predict across entire raster using focal windows.
        
        Parameters:
            raster_stack: Input raster stack
            batch_size: Number of windows to process in parallel
            show_progress: Whether to show progress bar
        
        Returns:
            Prediction raster with same spatial extent as input
        """
        height, width = raster_stack.data.shape[1:3]
        pad = self.window_size // 2
        
        # Pad raster to handle edges
        padded_data = jnp.pad(
            raster_stack.data,
            ((0, 0), (pad, pad), (pad, pad)),
            mode='reflect'
        )
        
        # Generate all window center coordinates
        rows, cols = jnp.meshgrid(
            jnp.arange(height),
            jnp.arange(width),
            indexing='ij'
        )
        coords = jnp.stack([rows.ravel(), cols.ravel()], axis=1)
        
        # Batch prediction
        n_pixels = coords.shape[0]
        predictions = []
        
        iterator = range(0, n_pixels, batch_size)
        if show_progress:
            iterator = tqdm(iterator, desc="Predicting")
        
        for start_idx in iterator:
            end_idx = min(start_idx + batch_size, n_pixels)
            batch_coords = coords[start_idx:end_idx]
            
            batch_preds = self._predict_batch(padded_data, batch_coords, pad)
            predictions.append(batch_preds)
        
        predictions = jnp.concatenate(predictions)
        return predictions.reshape(height, width)
    
    @partial(jit, static_argnums=(0,))
    def _predict_batch(
        self,
        padded_data: Float[Array, "bands h w"],
        coords: Float[Array, "batch 2"],
        pad: int
    ) -> Float[Array, "batch"]:
        """JIT-compiled batch prediction."""
        
        def predict_single(coord):
            r, c = coord
            # Extract window (accounting for padding offset)
            window = lax.dynamic_slice(
                padded_data,
                (0, r, c),
                (padded_data.shape[0], self.window_size, self.window_size)
            )
            # Reshape to (n_samples, n_features)
            # window is (bands, window_size, window_size)
            window_samples = window.reshape(window.shape[0], -1).T
            
            return self.predict_window(window_samples)
        
        return vmap(predict_single)(coords)


def predict_raster_parallel(
    predictor: FocalPredictor,
    raster_stack: RasterStack,
    n_blocks: int = 4
) -> Float[Array, "height width"]:
    """
    Parallel prediction using pmap across multiple devices.
    
    Splits the raster into blocks and processes in parallel.
    Handles edge collars to avoid boundary artifacts.
    """
    # Implementation for multi-GPU/TPU scenarios
    ...
```

### 5. High-Level API (`klrfome/__init__.py`)

```python
"""
KLRfome - Kernel Logistic Regression on Focal Mean Embeddings

A package for geographic distribution regression using kernel methods.
"""

from dataclasses import dataclass
from typing import Optional, Union, List
import jax.numpy as jnp

from .data.formats import TrainingData, RasterStack, SampleCollection
from .kernels.rbf import RBFKernel
from .kernels.rff import RandomFourierFeatures
from .kernels.distribution import MeanEmbeddingKernel
from .models.klr import KernelLogisticRegression, KLRFitResult
from .prediction.focal import FocalPredictor


@dataclass
class KLRfome:
    """
    High-level interface for KLRfome modeling.
    
    Example usage:
    
        # Initialize model
        model = KLRfome(
            sigma=1.0,
            lambda_reg=0.1,
            n_rff_features=256,
            window_size=5
        )
        
        # Prepare training data
        training_data = model.prepare_data(
            raster_stack=my_rasters,
            sites=site_geodataframe,
            n_background=1000,
            samples_per_location=20
        )
        
        # Fit model
        model.fit(training_data)
        
        # Predict
        prediction_raster = model.predict(prediction_rasters)
        
        # Save
        model.save_prediction("output.tif", prediction_raster)
    
    Parameters:
        sigma: RBF kernel bandwidth
        lambda_reg: KLR regularization strength
        n_rff_features: Number of random Fourier features (0 for exact kernel)
        window_size: Focal window size for prediction
        seed: Random seed for reproducibility
    """
    sigma: float = 1.0
    lambda_reg: float = 0.1
    n_rff_features: int = 256
    window_size: int = 3
    seed: int = 42
    
    # Fitted attributes (set after fit())
    _training_data: Optional[TrainingData] = None
    _similarity_matrix: Optional[jnp.ndarray] = None
    _fit_result: Optional[KLRFitResult] = None
    _distribution_kernel: Optional[MeanEmbeddingKernel] = None
    
    def __post_init__(self):
        # Initialize kernel
        if self.n_rff_features > 0:
            base_kernel = RandomFourierFeatures(
                sigma=self.sigma,
                n_features=self.n_rff_features,
                seed=self.seed
            )
        else:
            base_kernel = RBFKernel(sigma=self.sigma)
        
        self._distribution_kernel = MeanEmbeddingKernel(base_kernel)
        self._klr = KernelLogisticRegression(lambda_reg=self.lambda_reg)
    
    def prepare_data(
        self,
        raster_stack: Union[RasterStack, List[str]],
        sites: 'geopandas.GeoDataFrame',
        n_background: int = 1000,
        samples_per_location: int = 20,
        background_exclusion_buffer: Optional[float] = None,
        site_buffer: Optional[float] = None
    ) -> TrainingData:
        """
        Prepare training data from rasters and site locations.
        
        Parameters:
            raster_stack: RasterStack object or list of raster file paths
            sites: GeoDataFrame with site geometries
            n_background: Number of background sample locations
            samples_per_location: Samples to extract per site/background
            background_exclusion_buffer: Buffer around sites to exclude from background
            site_buffer: Buffer to apply around site points (if points, not polygons)
        
        Returns:
            TrainingData object ready for fitting
        """
        if isinstance(raster_stack, list):
            raster_stack = RasterStack.from_files(raster_stack)
        
        # Extract site samples
        site_collections = raster_stack.extract_at_points(
            sites,
            buffer_radius=site_buffer,
            n_samples=samples_per_location
        )
        for coll in site_collections:
            coll.label = 1
        
        # Generate background samples
        background_points = self._generate_background_points(
            raster_stack, sites, n_background, background_exclusion_buffer
        )
        background_collections = raster_stack.extract_at_points(
            background_points,
            n_samples=samples_per_location
        )
        for coll in background_collections:
            coll.label = 0
        
        return TrainingData(
            collections=site_collections + background_collections,
            feature_names=raster_stack.band_names,
            crs=raster_stack.crs
        )
    
    def fit(self, training_data: TrainingData) -> 'KLRfome':
        """
        Fit the KLRfome model.
        
        Parameters:
            training_data: Prepared TrainingData object
        
        Returns:
            self (for method chaining)
        """
        self._training_data = training_data
        
        # Build similarity matrix
        self._similarity_matrix = self._distribution_kernel.build_similarity_matrix(
            training_data.collections
        )
        
        # Fit KLR
        self._fit_result = self._klr.fit(
            self._similarity_matrix,
            training_data.labels
        )
        
        if not self._fit_result.converged:
            import warnings
            warnings.warn("KLR fitting did not converge")
        
        return self
    
    def predict(
        self,
        raster_stack: Union[RasterStack, List[str]],
        batch_size: int = 1000,
        show_progress: bool = True
    ) -> jnp.ndarray:
        """
        Generate predictions across a raster extent.
        
        Parameters:
            raster_stack: Rasters to predict on (must have same bands as training)
            batch_size: Windows to process per batch
            show_progress: Show progress bar
        
        Returns:
            2D array of predicted probabilities
        """
        if self._fit_result is None:
            raise RuntimeError("Model must be fit before prediction")
        
        if isinstance(raster_stack, list):
            raster_stack = RasterStack.from_files(raster_stack)
        
        predictor = FocalPredictor(
            distribution_kernel=self._distribution_kernel,
            klr_alpha=self._fit_result.alpha,
            training_data=self._training_data,
            window_size=self.window_size
        )
        
        return predictor.predict_raster(
            raster_stack,
            batch_size=batch_size,
            show_progress=show_progress
        )
    
    def save_prediction(
        self,
        path: str,
        prediction: jnp.ndarray,
        reference_raster: Optional[RasterStack] = None
    ):
        """Save prediction array as a GeoTIFF."""
        import rasterio
        
        if reference_raster is None:
            reference_raster = self._training_data  # Get from training
        
        # Implementation to write georeferenced output
        ...
    
    def cross_validate(
        self,
        training_data: TrainingData,
        n_folds: int = 5,
        stratified: bool = True
    ) -> dict:
        """
        Perform k-fold cross-validation.
        
        Returns:
            Dictionary with metrics per fold and aggregated statistics
        """
        ...
    
    def _generate_background_points(self, ...):
        """Generate random background sample locations."""
        ...


# Convenience aliases
__all__ = [
    'KLRfome',
    'TrainingData',
    'RasterStack', 
    'SampleCollection',
    'RBFKernel',
    'RandomFourierFeatures',
    'MeanEmbeddingKernel',
    'KernelLogisticRegression',
    'FocalPredictor',
]
```

---

## Testing Requirements

### Unit Tests (`tests/`)

```
tests/
├── conftest.py              # Shared fixtures
├── test_kernels.py          # Kernel computations
├── test_klr.py              # KLR fitting
├── test_prediction.py       # Focal prediction
├── test_data.py             # Data structures
└── test_integration.py      # End-to-end tests
```

#### Key Test Cases

```python
# tests/test_kernels.py

def test_rff_approximates_rbf():
    """Random Fourier Features should approximate exact RBF."""
    X = jnp.array([[0, 0], [1, 0], [0, 1], [1, 1]])
    sigma = 1.0
    
    exact = RBFKernel(sigma=sigma)(X, X)
    approx = RandomFourierFeatures(sigma=sigma, n_features=1000)(X, X)
    
    assert jnp.allclose(exact, approx, atol=0.1)


def test_rbf_kernel_properties():
    """RBF kernel should satisfy kernel properties."""
    X = random.normal(key, (10, 5))
    K = RBFKernel(sigma=1.0)(X, X)
    
    # Symmetric
    assert jnp.allclose(K, K.T)
    
    # Positive semi-definite
    eigenvalues = jnp.linalg.eigvalsh(K)
    assert jnp.all(eigenvalues >= -1e-10)
    
    # Diagonal is 1
    assert jnp.allclose(jnp.diag(K), 1.0)


def test_distribution_kernel_symmetry():
    """Distribution kernel should be symmetric."""
    samples_a = random.normal(key1, (20, 5))
    samples_b = random.normal(key2, (15, 5))
    
    kernel = MeanEmbeddingKernel(RBFKernel(sigma=1.0))
    
    assert jnp.isclose(kernel(samples_a, samples_b), kernel(samples_b, samples_a))


# tests/test_klr.py

def test_klr_fits_separable_data():
    """KLR should achieve high accuracy on separable data."""
    # Generate clearly separable distributions
    ...
    
    result = klr.fit(K, y)
    probs = klr.predict_proba(K, result.alpha)
    accuracy = jnp.mean((probs > 0.5) == y)
    
    assert accuracy > 0.9
    assert result.converged


def test_klr_regularization_effect():
    """Higher lambda should produce smaller alpha magnitudes."""
    klr_low = KernelLogisticRegression(lambda_reg=0.01)
    klr_high = KernelLogisticRegression(lambda_reg=10.0)
    
    result_low = klr_low.fit(K, y)
    result_high = klr_high.fit(K, y)
    
    assert jnp.linalg.norm(result_high.alpha) < jnp.linalg.norm(result_low.alpha)


# tests/test_integration.py

def test_end_to_end_simulated():
    """Full pipeline should work on simulated data."""
    # Generate simulated rasters and sites
    ...
    
    model = KLRfome(sigma=1.0, lambda_reg=0.1, n_rff_features=256)
    training_data = model.prepare_data(rasters, sites, n_background=100)
    model.fit(training_data)
    predictions = model.predict(rasters)
    
    assert predictions.shape == (rasters.height, rasters.width)
    assert jnp.all((predictions >= 0) & (predictions <= 1))
```

---

## Documentation Requirements

1. **README.md**: Installation, quickstart, basic usage
2. **API Reference**: Auto-generated from docstrings (sphinx/mkdocs)
3. **Tutorial Notebooks**:
   - `01_quickstart.ipynb`: Basic usage with simulated data
   - `02_real_data.ipynb`: Working with real rasters and shapefiles
   - `03_hyperparameter_tuning.ipynb`: Cross-validation and tuning
   - `04_advanced_kernels.ipynb`: Customizing kernels
4. **Theory Guide**: Mathematical background (can adapt from original JOSS paper)

---

## Implementation Priorities

### Phase 1: Core Functionality (MVP)
1. Data structures (`SampleCollection`, `TrainingData`, `RasterStack`)
2. RBF kernel (exact)
3. Random Fourier Features
4. Mean embedding distribution kernel
5. KLR with IRLS
6. Basic focal prediction (single-threaded)
7. High-level `KLRfome` class
8. Basic tests

### Phase 2: Performance
1. JIT compilation for all hot paths
2. Batched focal prediction with `vmap`
3. GPU support verification
4. Benchmarking against R implementation

### Phase 3: Usability
1. Rasterio integration for I/O
2. GeoPandas integration for sites
3. Progress bars
4. Visualization utilities
5. Cross-validation
6. Documentation and tutorials

### Phase 4: Extensions
1. Additional kernel options
2. Multi-scale focal windows
3. Parallel prediction across devices
4. Model serialization/persistence

---

## Notes for AI Assistant

1. **Prioritize correctness over cleverness**: The IRLS algorithm and kernel computations must be numerically stable. Use established formulations.

2. **JAX idioms**: Use `jit`, `vmap`, and functional style throughout. Avoid Python loops where vectorization is possible.

3. **Type hints**: Use `jaxtyping` for array shapes. This aids both documentation and catches shape errors early.

4. **Testing**: Property-based testing with Hypothesis is valuable for numerical code. Test that kernels satisfy mathematical properties (symmetry, PSD, etc.).

5. **Memory efficiency**: The full kernel matrix can be large. The RFF approach avoids this, but when using exact kernels, be mindful of memory.

6. **Reproducibility**: All random operations should accept a seed/key. Default to deterministic behavior.

7. **Geospatial conventions**: Follow rasterio conventions for transforms, CRS handling, etc. Don't reinvent the wheel for I/O.

8. **Error messages**: Provide helpful error messages when shapes don't match, data is invalid, etc. Users of this package may not be ML experts.
