"""
Setup script for spaGAPA package.
For modern installations, use pyproject.toml instead.
"""

from setuptools import setup, find_packages

# Read version from __init__.py
with open("spagapa/__init__.py", "r") as f:
    for line in f:
        if line.startswith("__version__"):
            version = line.split("=")[1].strip().strip('"').strip("'")
            break

# Read long description from README
try:
    with open("README.md", "r", encoding="utf-8") as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = "spaGAPA: spatial Gaussian process-based APA analyzer"

setup(
    name="spagapa",
    version=version,
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "pandas>=1.3.0",
        "scikit-learn>=1.0.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "plotly>=5.0.0",
        "pysam>=0.19.0",
        "anndata>=0.8.0",
        "scanpy>=1.9.0",
        "joblib>=1.1.0",
        "tqdm>=4.62.0",
        "statsmodels>=0.13.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "hypothesis>=6.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
            "mypy>=0.950",
            "isort>=5.10.0",
        ],
        "docs": [
            "sphinx>=4.5.0",
            "sphinx-rtd-theme>=1.0.0",
            "nbsphinx>=0.8.0",
            "jupyter>=1.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "spagapa=spagapa.cli:main",
        ],
    },
    python_requires=">=3.8",
    author="spaGAPA Development Team",
    author_email="your.email@example.com",
    description="spatial Gaussian process-based APA analyzer",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/spaGAPA",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
