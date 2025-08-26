import torch
import torch.nn.functional as F
import numpy as np
from typing import Tuple
import pytest
# from maldi_nn.utils import topf as _topf


class PadSequence(object):
    """ Pad the input sequence to the target shape with the padding value. All 
    padding is on the end of the tensor.

    Args:
        shape (Tuple[int,]): The target shape of the output tensor.
        padding_value (int): The value for the padded elements.

    Returns:
        Tensor: The output tensor with the target shape.

    Example:
        >>> transformer = PadSequence((5,))
        >>> input_tensor = torch.tensor([1, 2, 3])
        >>> transformer(input_tensor)
        tensor([1, 2, 3, 0, 0]) 
    """
    def __init__(self, shape:Tuple[int,], padding_value=0):
        self.shape = shape
        self.padding_value = padding_value

    def __call__(self, vector):
        # Create an output tensor filled with the padding value, with the target shape
        output_vector = torch.full(self.shape, self.padding_value)
        
        # Create slices for each dimension based on the smaller of the output shape and input vector's shape
        slices = tuple(slice(0, min(v, o)) for v, o in zip(vector.shape, self.shape))
        
        # Assign the input vector into the corresponding slice of the output vector
        output_vector[slices] = vector
        
        return output_vector
    
class BinSpectrum(object):
    """Bin the input spectrum into m/z bins of fixed width using PyTorch. 
    Intensity values within each bin are summed.

    Args:
        bin_width (float): The width of each bin.
        min_mz (float): The minimum m/z value.
        max_mz (float): The maximum m/z value.

    Returns:
        Tensor: The binned spectrum.
    """

    def __init__(self, bin_width: float, min_mz: float, max_mz: float):
        self.bin_width = bin_width
        self.min_mz = min_mz
        self.max_mz = max_mz
        self.n_bins = int((self.max_mz - self.min_mz) / self.bin_width)
    
    def __call__(self, spectrum: torch.Tensor) -> torch.Tensor:
        # spectrum: Tensor of shape (N, 2) with m/z and intensity
        mz = spectrum[:, 0]
        intensity = spectrum[:, 1]

        # Compute bin indices for each m/z value
        bin_indices = ((mz - self.min_mz) / self.bin_width).floor().long()

        # Mask out-of-range indices
        valid_mask = (bin_indices >= 0) & (bin_indices < self.n_bins)
        bin_indices = bin_indices[valid_mask]
        intensity = intensity[valid_mask]

        # Sum intensities within bins
        binned_spectrum = torch.zeros(self.n_bins, dtype=spectrum.dtype, device=spectrum.device)
        binned_spectrum.index_add_(0, bin_indices, intensity)

        return binned_spectrum
    
class SquareRootTransform(object):
    """Performs a square root transformation on the intensities of the spectrum.
    
    Args:
        None
        
    Returns:
        torch.Tensor: The transformed spectrum.
    """
    def __init__(self):
        pass

    def __call__(self, spectrum):
        if spectrum.ndim == 1:
            spectrum[1] = torch.sqrt(spectrum[1])
        elif spectrum.ndim == 2:
            spectrum[:, 1] = torch.sqrt(spectrum[:, 1])
        else:
            raise ValueError(f"Expected a 1D or 2D tensor with m/z and intensity values. Instead got {spectrum.shape}")
        return spectrum

class NormalizeIntensity(object):
    """Normalizes intensity values along a specified dimension using the Euclidean norm,
    while preserving other dimensions such as m/z.
    
    Args:
        dim (int): The dimension to normalize (default is -1, the last dimension).
        
    Returns:
        torch.Tensor: The tensor with normalized intensities along the specified dimension.
    """
    
    def __init__(self):
        pass

    def __call__(self, spectrum: torch.Tensor) -> torch.Tensor:
        if spectrum.ndim == 1:
            if torch.sum(spectrum) == 0:
                return spectrum

            norm = torch.norm(spectrum, p=2)
            spectrum = spectrum / norm
            return spectrum

        if spectrum.ndim == 2:
            if torch.sum(spectrum[:, 1]) == 0:
                return spectrum
                raise ValueError(f"Expected a 2D tensor with at least two columns (m/z and intensity). Got shape {spectrum.shape}")

            if spectrum.shape[1] != 2:
                raise ValueError(f"Expected a 2D tensor with m/z and intensity values. Got shape {spectrum.shape}")

            norm = torch.norm(spectrum[:, 1], p=2)
            spectrum[:, 1] = spectrum[:, 1] / norm
            return spectrum

        raise ValueError(f"Expected a 1D or 2D tensor with m/z and intensity values. Got shape {spectrum.shape}")

class L1NormalizeIntensity(object):
    """Normalizes intensity values along a specified dimension using the L1 norm,
    while preserving other dimensions such as m/z.
    
    Args:
        dim (int): The dimension to normalize (default is -1, the last dimension).
        
    Returns:
        torch.Tensor: The tensor with normalized intensities along the specified dimension.
    """
    
    def __init__(self):
        pass

    def __call__(self, spectrum: torch.Tensor) -> torch.Tensor:
        if spectrum.ndim == 1:
            if torch.sum(spectrum) == 0:
                return spectrum

            norm = torch.norm(spectrum, p=1)
            spectrum = spectrum / norm
            return spectrum

        if spectrum.ndim == 2:
            if torch.sum(spectrum[:, 1]) == 0:
                return spectrum
                raise ValueError(f"Expected a 2D tensor with at least two columns (m/z and intensity). Got shape {spectrum.shape}")

            if spectrum.shape[1] != 2:
                raise ValueError(f"Expected a 2D tensor with m/z and intensity values. Got shape {spectrum.shape}")

            norm = torch.norm(spectrum[:, 1], p=1)
            spectrum[:, 1] = spectrum[:, 1] / norm
            return spectrum

        raise ValueError(f"Expected a 1D or 2D tensor with m/z and intensity values. Got shape {spectrum.shape}")   

class BinarizeIntensity(object):
    """Binarizes the intensity values of the spectrum. If the intensity is greater than threshold, it is set to 1; otherwise, set it to 0.
    
    Args:
        threshold (float): The threshold for binarization.

    Returns:
        torch.Tensor: The binarized spectrum.
    """
    def __init__(self, threshold: float = 0.02):
        self.threshold = threshold

    def __call__(self, spectrum):
        if len(spectrum.shape) == 1:
            spectrum[:] = (spectrum > self.threshold).float()
        elif len(spectrum.shape) == 2:
            spectrum[:, 1] = (spectrum[:, 1] > self.threshold).float()
        else:
            raise ValueError(f"Expected a 1D or 2D array with m/z and intensity values. Instead got {spectrum.shape}")
        return spectrum

class ScaleIntensity(object):
    """Ensures the max intensity of the spectrum is 1.0. 

    Args:
        None

    Returns:
        torch.tensor: The scaled spectrum.
    """
    def __init__(self):
        pass

    def __call__(self, spectrum):
        
        if len(spectrum.shape) == 1:
            if sum(spectrum) == 0:
                raise ValueError("Empty spectrum")
            spectrum[1] = spectrum[1] / torch.max(spectrum[1])
        elif len(spectrum.shape) == 2:
            if spectrum.shape[1] < 1:
                raise ValueError(f"Expected a 2D array with at least two columns (m/z and intensity). Instead got {spectrum.shape}")            
            spectrum[:, 1] = spectrum[:, 1] / torch.max(spectrum[:, 1])
        else:
            raise ValueError(f"Expected a 1D or 2D array with m/z and intensity values. Instead got {spectrum.shape}")
        return spectrum
    
class SelectMassRange(object):
    """ Selects the mass range of a spectrum to the specified range. Both endpoints are inclusive.

    Args:
        min_mz (float): The minimum m/z value.
        max_mz (float): The maximum m/z value.

    Returns:
        np.ndarray: The reduced spectrum.
    """

    def __init__(self, min_mz:float, max_mz:float):
        self.min_mz = min_mz
        self.max_mz = max_mz

    def __call__(self, spectrum):
        if len(spectrum.shape) != 2:
            raise ValueError(f"Expected a 2D array with m/z and intensity values. Instead got {spectrum.shape}")
        
        # Filter the spectrum based on the m/z values
        mask = (spectrum[:, 0] >= self.min_mz) & (spectrum[:, 0] <= self.max_mz)
        spectrum = spectrum[mask]
        return spectrum
    
# class SelectTopKPeaks(object):
#     """Select the top peaks by intensity from a spectrum.
    
#     Args:
#         k (int): The number of peaks to select.
        
#     Returns:
#         torch.Tensor: The reduced spectrum.
#     """
#     def __init__(self, k:int):
#         self.k = k

#     def __call__(self, spectrum):
#         if len(spectrum.shape) != 2:
#             raise ValueError(f"Expected a 2D tensor with m/z and intensity values. Instead got {spectrum.shape}")
        
#         if spectrum.shape[1] < 2:
#             raise ValueError(f"Expected a 2D tensor with at least two columns (m/z and intensity). Instead got {spectrum.shape}")
        
#         spectrum = np.array(spectrum)

#         # Sort the spectrum by intensity (second column) and get indices
#         sorted_indices = np.argsort(spectrum[:, 1])
        
#         # Select the top k peaks by intensity
#         spectrum = spectrum[sorted_indices[-self.k:]]

#         # Sort by m/z values
#         sorted_indices = np.argsort(spectrum[:, 0])
#         spectrum = spectrum[sorted_indices]
        
#         return spectrum
    
class SelectTopKPeaks:
    def __init__(self, k: int):
        self.k = k

    def __call__(self, spectrum: torch.Tensor) -> torch.Tensor:
        if spectrum.ndim != 2 or spectrum.shape[1] < 2:
            raise ValueError(f"Expected a 2D tensor with at least two columns (m/z and intensity). Got shape {spectrum.shape}")

        k = min(self.k, spectrum.shape[0])  # handle case where fewer peaks than k

        # Select top-k by intensity (second column)
        intensities = spectrum[:, 1]
        topk_values, topk_indices = torch.topk(intensities, k=k, largest=True, sorted=False)
        topk_peaks = spectrum[topk_indices]

        # Sort selected peaks by m/z (first column)
        sorted_indices = torch.argsort(topk_peaks[:, 0])
        sorted_topk_peaks = topk_peaks[sorted_indices]

        return sorted_topk_peaks
    
class ExcludeMassRange(object):
    """ Excludes the mass range of a spectrum to the specified range. Both endpoints are inclusive.
    
    Args:
        min_mz (float): The minimum m/z value.
        max_mz (float): The maximum m/z value.
        
    Returns:
        np.ndarray: The reduced spectrum.
    """
    def __init__(self, min_mz:float, max_mz:float):
        self.min_mz = min_mz
        self.max_mz = max_mz

    def __call__(self, spectrum):
        if len(spectrum.shape) != 2:
            raise ValueError(f"Expected a 2D array with m/z and intensity values. Instead got {spectrum.shape}")
        
        # Filter the spectrum based on the m/z values
        mask = (spectrum[:, 0] < self.min_mz) | (spectrum[:, 0] > self.max_mz)
        spectrum = spectrum[mask]
        return spectrum
    
class PadToLength:
    """Pads a tensor along a specified dimension to a fixed length.
    
    Args:
        length (int): The target length of the specified dimension.
        dim (int): The dimension to pad (default is 0).
        padding_value (float): The value for the padded elements (default is -1.0).
        
    Returns:
        torch.Tensor: The padded (or truncated) tensor.
    """
    def __init__(self, length: int, dim: int = 0, padding_value: float = -1.0):
        self.length = length
        self.dim = dim
        self.padding_value = padding_value

    def __call__(self, sequence: torch.Tensor) -> torch.Tensor:
        dim = self.dim if self.dim >= 0 else sequence.ndim + self.dim
        current_length = sequence.shape[dim]

        if current_length >= self.length:
            return torch.narrow(sequence, dim, 0, self.length)

        pad_size = self.length - current_length
        pad_shape = list(sequence.shape)
        pad_shape[dim] = pad_size
        padding = torch.full(pad_shape, self.padding_value, dtype=sequence.dtype, device=sequence.device)

        return torch.cat([sequence, padding], dim=dim)
    
class NoiseInjection(object):
    """
    Args:
        noise_factor (float): Max value of the noise to be added to the spectrum.
    """
    def __init__(self, noise_factor=0.05):
        self.noise_factor = noise_factor

    def __call__(self, spectrum):
        if len(spectrum.shape) == 1:
            spectrum = np.array(spectrum)
            noise = np.random.uniform(0, self.noise_factor, spectrum.shape)
            spectrum += noise
            return spectrum
        elif len(spectrum.shape) == 2:
            max_mz = spectrum[:, 0].max()
            min_mz = spectrum[:, 0].min()

            noise_quantity = 0.25 * sum(spectrum[:, 1] > self.noise_factor)

            noise_mz = np.random.uniform(min_mz, max_mz, int(noise_quantity))
            noise_intensity = np.random.uniform(0, self.noise_factor, int(noise_quantity))
            noise = np.column_stack((noise_mz, noise_intensity))
            spectrum = np.vstack((spectrum, noise))
            # Sort by m/z values (may matter for some models)
            sorted_indices = np.argsort(spectrum[:, 0])
            spectrum = spectrum[sorted_indices]
            return spectrum
            

# class topf(object):
#     """Applies the topf transformation described in TODO (see package documentation).
    
#     Args:
#         None

#     Returns:
#         np.ndarray: The transformed spectrum.
#     """
#     def __init__(self):
#         pass

#     def __call__(self, spectrum):
#         spectrum = np.array(spectrum)

#         assert len(spectrum.shape) == 2, f"Expected a 2D array with m/z and intensity values. Instead got {spectrum.shape}"
#         assert spectrum.shape[1] == 2, f"Expected a 2D array with two rows (m/z and intensity). Instead got {spectrum.shape}."

#         persistance = _topf.PersistenceTransformer().fit_transform(spectrum)

#         return persistance

class toTensor(object):
    """Converts a numpy array to a PyTorch tensor.
    
    Args:
        None

    Returns:
        torch.Tensor: The PyTorch tensor.
    """
    def __init__(self):
        pass

    def __call__(self, array):
        return torch.tensor(array).to(torch.float32)

@pytest.fixture
def test_spectrum():
    yield np.array([[1, 20], [2, 30], [3, 40], [4, 50], [5, 60]])   # m/z, intensity

def test_bin_spectrum(test_spectrum):
    transformer = BinSpectrum(3, 1, 6)
    output = transformer(test_spectrum)

    assert output.shape == (2,)
    assert output[0] == 20 + 30 + 40, f"First bin should contain the sum of intensities 20, 30, and 40, but got {output[0]}"
    assert output[1] == 50 + 60, f"Second bin should contain the sum of intensities 50 and 60, but got {output[1]}"

def test_reduce_mass_range(test_spectrum):
    transformer = SelectMassRange(2, 4)
    output = transformer(test_spectrum)

    assert output.shape == (3, 2), f"Expected shape (3, 2), but got {output.shape}"
    assert np.all(output[:, 0] >= 2), "All m/z values should be greater than or equal to 2"
    assert np.all(output[:, 0] <= 4), "All m/z values should be less than or equal to 4"

def test_exclude_mass_range(test_spectrum):
    transformer = ExcludeMassRange(2, 4)
    output = transformer(test_spectrum)

    assert output.shape == (2, 2), f"Expected shape (2, 2), but got {output.shape}"
    assert np.all((output[:, 0] < 2) | (output[:, 0] > 4)), "All m/z values should be less than 2 or greater than 4"