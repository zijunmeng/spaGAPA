Analysis presets
================

Four user-facing modes (:mod:`spagapa.presets`):

===================  =========================================================
Preset               When to use
===================  =========================================================
``auto``             Default. Profiles the matrix (n_spots, observed
                    fraction, zero-encoding) and resolves to ``standard``
                    or ``highres_fast`` automatically.
``standard``         Visium-scale: thousands of spots, moderate sparsity.
                    Dense/exact GP defaults.
``highres_accuracy`` Large chips where fidelity matters more than wall time.
``highres_fast``     Large chips, production runs. Sparse GP + chunked
                    factorizer + Leiden domains. The mode used for all
                    Stereo-seq results in the manuscript.
===================  =========================================================

The resolver treats data as high-resolution when any of:

* ``n_spots >= 1000`` and observed fraction < 0.35, or
* zero-encoded APA-index matrices with sparse positive signal
  (details in ``presets.profile_spatial_apa_matrix``).

What each preset switches
-------------------------

=========================  ==================  ==================
Capability                 standard            highres_*
=========================  ==================  ==================
GP backend                 exact GP            sparse GP (m inducing)
Inducing points            n/a                 ``n_inducing`` (200 default)
BioML graph factorizer     in-memory           chunked (n_spots > 20k)
Domain detection           spectral/kmeans     Leiden (n_spots > 5k)
=========================  ==================  ==================

Overriding preset internals
---------------------------

Presets set defaults; every knob remains a parameter (``SpaGAPA(...)`` or
CLI flags). Common high-resolution tuning:

.. code-block:: python

    SpaGAPA(
        analysis_preset="highres_fast",
        n_inducing=200,
        sparse_gp_length_scale="auto",   # coordinate-adaptive
        sparse_gp_noise_level=0.1,
        bioml_domain_method="spectral",
        highres_bioml_gp_blend=0.1,
        highres_bioml_spatial_weight=0.1,
        highres_bioml_expression_weight=0.7,
        highres_bioml_apa_weight=0.2,
    )

The resolved preset actually used is reported in the run summary
(``results['analysis_preset']``) — check it when passing ``auto``.
