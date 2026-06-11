"""
Unit tests for scAPAtrap wrapper.
"""

import pytest
import subprocess
from unittest.mock import Mock, patch, MagicMock
from spagapa.io.scapatrap_wrapper import ScAPAtrapWrapper


class TestScAPAtrapWrapper:
    """Test suite for ScAPAtrapWrapper class."""
    
    def test_init(self):
        """Test wrapper initialization."""
        with patch.object(ScAPAtrapWrapper, '_check_r_installation'):
            with patch.object(ScAPAtrapWrapper, '_check_scapatrap_installation'):
                wrapper = ScAPAtrapWrapper()
                assert wrapper.r_executable == "Rscript"
    
    def test_check_r_installation_success(self):
        """Test R installation check when R is available."""
        with patch('subprocess.run') as mock_run:
            mock_run.return_value = Mock(
                stdout="R version 4.0.0",
                returncode=0
            )
            with patch.object(ScAPAtrapWrapper, '_check_scapatrap_installation'):
                wrapper = ScAPAtrapWrapper()
                # Should not raise exception
    
    def test_check_r_installation_failure(self):
        """Test R installation check when R is not available."""
        with patch('subprocess.run') as mock_run:
            mock_run.side_effect = FileNotFoundError()
            with pytest.raises(RuntimeError, match="R executable not found"):
                ScAPAtrapWrapper()
    
    def test_generate_r_script(self):
        """Test R script generation."""
        with patch.object(ScAPAtrapWrapper, '_check_r_installation'):
            with patch.object(ScAPAtrapWrapper, '_check_scapatrap_installation'):
                wrapper = ScAPAtrapWrapper()
                
                script = wrapper._generate_r_script(
                    bam_file="test.bam",
                    output_dir="output",
                    genome_fasta="genome.fa",
                    gtf_file="genes.gtf",
                    barcode_file=None,
                    tails_search="genome",
                    n_cores=4
                )
                
                assert "library(scAPAtrap)" in script
                assert "test.bam" in script
                assert "genome.fa" in script
                assert "genes.gtf" in script


@pytest.mark.integration
class TestScAPAtrapIntegration:
    """Integration tests for scAPAtrap wrapper (requires R and scAPAtrap)."""
    
    def test_r_available(self):
        """Test if R is available in the system."""
        try:
            result = subprocess.run(
                ["Rscript", "--version"],
                capture_output=True,
                text=True,
                check=True
            )
            assert "R scripting front-end" in result.stdout or "version" in result.stdout
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("R not available in system")
    
    def test_scapatrap_available(self):
        """Test if scAPAtrap is installed in R."""
        r_code = 'if (require("scAPAtrap", quietly = TRUE)) cat("OK")'
        try:
            result = subprocess.run(
                ["Rscript", "-e", r_code],
                capture_output=True,
                text=True,
                check=True
            )
            if "OK" not in result.stdout:
                pytest.skip("scAPAtrap not installed in R")
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("Cannot check scAPAtrap installation")
