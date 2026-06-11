"""
Example 09: Benchmark — spaGAPA vs baseline methods

Runs the full simulation benchmark and generates publication-quality figures.
Results are saved to benchmark_results/.
"""

import numpy as np
from spagapa.benchmark import (
    run_full_benchmark,
    run_mob_benchmark,
    generate_all_benchmark_figures,
)


def main():
    print("=" * 60)
    print("spaGAPA Benchmark")
    print("=" * 60)

    # ── 1. Simulation benchmark (3 dropout scenarios) ──────────────────────────
    print("\n[1/2] Simulation benchmark (Mean / Median / KNN / spaGAPA-GP)")
    results = run_full_benchmark(
        simulator_kwargs=dict(n_spots=200, n_genes=40, n_domains=4,
                              noise_level=0.08),
        scenarios=[
            {'name': 'low_dropout',    'dropout_rate': 0.20},
            {'name': 'medium_dropout', 'dropout_rate': 0.50},
            {'name': 'high_dropout',   'dropout_rate': 0.70},
        ],
        output_dir='benchmark_results',
        random_state=42
    )

    print("\nGenerating figures...")
    figs = generate_all_benchmark_figures(
        results,
        output_dir='benchmark_results',
        dropout_rates=[0.20, 0.50, 0.70]
    )
    print(f"Saved {len(figs)} figures to benchmark_results/")

    # ── 2. MOB benchmark (high-fidelity simulator) ─────────────────────────────
    print("\n[2/2] MOB benchmark (concentric-ring simulator)")
    mob_out = run_mob_benchmark(
        data_dir=None,          # set to real data path when available
        output_dir='benchmark_results/mob',
        n_spots=200,
        n_genes=40,
        random_state=42
    )

    print("\n" + "=" * 60)
    print("Benchmark complete.")
    print(f"  Simulation figures : benchmark_results/")
    print(f"  MOB figures        : benchmark_results/mob/")
    print("=" * 60)


if __name__ == '__main__':
    main()
