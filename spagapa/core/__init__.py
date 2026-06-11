"""
Core data structures for spaGAPA.

This module provides the fundamental data structures for APA analysis:
- APADataset: Main data container (wraps AnnData)
- APASite: Individual APA site representation
- APASiteCollection: Collection of APA sites with filtering
"""

from spagapa.core.apa_dataset import APADataset
from spagapa.core.apa_site import APASite, APASiteCollection

__all__ = [
    "APADataset",
    "APASite",
    "APASiteCollection",
]
