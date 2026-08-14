"""
compute surface pressure core module

This module contains the core functions for computing surface pressure from FWH data. It includes functions for extracting files, processing data, and generating results.

"""
from .extractor import extract_files, extract_surface, extract_data, extract_surface_line
from .surface_psd import PSD_surface_data, compute_PSD_block, source_psd
from .surface_csd import CSD_surface_data, compute_CSD_block, source_csd
from .utils import (
    timer,
    next_greater_power_of_2,
    fft_surface_data,
    spod_parser,
    hammwin,
    replace_zeros_vectorized,
)
from .ComputeSurfacePressure import SurfacePressure, parse_arguments, main
 
 
__all__ = [
    "extract_files",
    "extract_surface",
    "extract_data",
    "extract_surface_line",
    "PSD_surface_data",
    "compute_PSD_block",
    "source_psd",
    "CSD_surface_data",
    "compute_CSD_block",
    "source_csd",
    "timer",
    "next_greater_power_of_2",
    "fft_surface_data",
    "spod_parser",
    "hammwin",
    "replace_zeros_vectorized",
    "SurfacePressure",
    "parse_arguments",
    "main",
]
 
__version__ = "0.1.0"