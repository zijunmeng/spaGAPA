"""
Example: Spatial-aware APA calling

This example demonstrates how to use the spatial validation and quality
filtering modules to identify high-confidence APA sites from spatial
transcriptomics data.
"""

import numpy as np
from spagapa.core import APASite, APASiteCollection
from spagapa.spatial import SpatialNeighbors
from spagapa.calling import SpatialValidator, QualityFilter

# Simulate spatial coordinates (e.g., 10x Visium spots)
np.random.seed(42)
n_spots = 100
spatial_coords = np.random.rand(n_spots, 2) * 10  # 10x10 grid

# Simulate candidate APA sites
candidate_sites = APASiteCollection([
    APASite(
        chr='chr1',
        start=i*1000,
        end=i*1000+50,
        strand='+',
        gene_id=f'gene{i}',
        gene_name=f'GENE{i}'
    )
    for i in range(20)
])

# Simulate APA count matrix (20 genes x 100 spots)
# Some genes have spatially coherent patterns, others are noisy
apa_counts = np.zeros((20, n_spots))

# Gene 0-5: Strong spatial pattern (clustered in one region)
for i in range(6):
    # Spots in upper-left quadrant have high counts
    mask = (spatial_coords[:, 0] < 5) & (spatial_coords[:, 1] < 5)
    apa_counts[i, mask] = np.random.poisson(10, size=mask.sum())

# Gene 6-10: Moderate spatial pattern
for i in range(6, 11):
    # Spots in lower-right quadrant have moderate counts
    mask = (spatial_coords[:, 0] > 5) & (spatial_coords[:, 1] > 5)
    apa_counts[i, mask] = np.random.poisson(5, size=mask.sum())

# Gene 11-15: Weak/noisy pattern (random spots)
for i in range(11, 16):
    random_spots = np.random.choice(n_spots, size=10, replace=False)
    apa_counts[i, random_spots] = np.random.poisson(3, size=10)

# Gene 16-19: Very sparse (only 1-2 spots)
for i in range(16, 20):
    random_spots = np.random.choice(n_spots, size=2, replace=False)
    apa_counts[i, random_spots] = np.random.poisson(2, size=2)

print("=" * 60)
print("Spatial-Aware APA Calling Example")
print("=" * 60)
print(f"\nInput:")
print(f"  - {len(candidate_sites)} candidate APA sites")
print(f"  - {n_spots} spatial spots")
print(f"  - {np.sum(apa_counts > 0)} non-zero counts")

# Step 1: Spatial validation
print("\n" + "=" * 60)
print("Step 1: Spatial Validation")
print("=" * 60)

validator = SpatialValidator(
    n_neighbors=6,
    support_threshold=0.3,
    method='knn'
)
validator.fit(spatial_coords)

validated_sites, support_scores = validator.validate_sites(
    candidate_sites,
    apa_counts
)

print(f"\nSpatial validation results:")
print(f"  - Sites passing validation: {len(validated_sites)}/{len(candidate_sites)}")
print(f"  - Mean support score: {np.mean(support_scores):.3f}")
print(f"  - Support score range: [{np.min(support_scores):.3f}, {np.max(support_scores):.3f}]")

# Show which sites passed
print(f"\nSites with high spatial support (>0.5):")
for i, (site, score) in enumerate(zip(candidate_sites, support_scores)):
    if score > 0.5:
        print(f"  - {site.gene_name}: support={score:.3f}")

# Step 2: Quality filtering
print("\n" + "=" * 60)
print("Step 2: Quality Filtering")
print("=" * 60)

qf = QualityFilter(
    min_read_count=20,
    min_spots=5,
    min_mean_count=2.0,
    min_spatial_support=0.3
)

filtered_sites, pass_mask = qf.filter_sites(
    candidate_sites,
    apa_counts,
    spatial_support=support_scores
)

print(f"\nQuality filtering results:")
print(f"  - Sites passing all filters: {len(filtered_sites)}/{len(candidate_sites)}")
print(f"  - Filters applied:")
print(f"    * min_read_count: {qf.min_read_count}")
print(f"    * min_spots: {qf.min_spots}")
print(f"    * min_mean_count: {qf.min_mean_count}")
print(f"    * min_spatial_support: {qf.min_spatial_support}")

# Step 3: Generate QC report
print("\n" + "=" * 60)
print("Step 3: Quality Control Report")
print("=" * 60)

qc_report = qf.generate_qc_report(
    candidate_sites,
    apa_counts,
    spatial_support=support_scores
)

print(f"\nQC report summary:")
print(f"  - Total sites: {len(qc_report)}")
print(f"  - Sites passing all filters: {qc_report['pass_all'].sum()}")
print(f"\nFilter pass rates:")
print(f"  - Read count: {qc_report['pass_read_count'].sum()}/{len(qc_report)}")
print(f"  - Spot count: {qc_report['pass_spot_count'].sum()}/{len(qc_report)}")
print(f"  - Mean count: {qc_report['pass_mean_count'].sum()}/{len(qc_report)}")
print(f"  - Spatial support: {qc_report['pass_spatial_support'].sum()}/{len(qc_report)}")

# Show top 5 sites by quality
print(f"\nTop 5 sites by total count:")
top_sites = qc_report.nlargest(5, 'total_count')[
    ['gene_id', 'total_count', 'n_spots', 'mean_count', 'spatial_support', 'pass_all']
]
print(top_sites.to_string(index=False))

# Step 4: Spatial consistency filtering
print("\n" + "=" * 60)
print("Step 4: Spatial Consistency Filtering")
print("=" * 60)

filtered_counts = validator.filter_by_spatial_consistency(
    apa_counts,
    min_support=0.3
)

n_filtered = np.sum((apa_counts > 0) & (filtered_counts == 0))
total_nonzero = np.sum(apa_counts > 0)

print(f"\nSpatial consistency filtering:")
print(f"  - Original non-zero counts: {total_nonzero}")
print(f"  - Filtered counts: {n_filtered}")
print(f"  - Remaining counts: {np.sum(filtered_counts > 0)}")
print(f"  - Filtering rate: {100*n_filtered/total_nonzero:.1f}%")

# Summary
print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"\nFinal high-confidence APA sites: {len(filtered_sites)}")
print(f"  - Spatially validated: ✓")
print(f"  - Quality filtered: ✓")
print(f"  - Ready for downstream analysis")

print("\n" + "=" * 60)
print("Example completed successfully!")
print("=" * 60)
