"""Check repeated timings and plot stable/full runs apart from the aborted witness."""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/funtides-lts-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def analyze(root):
    protocol = json.loads((root / 'timing_protocol.json').read_text())
    result = json.loads((root / 'test_result.json').read_text())
    raw = np.genfromtxt(root / 'timing_samples.csv', delimiter=',', names=True, dtype=None, encoding='utf8')
    summary = np.genfromtxt(root / 'timing_summary.csv', delimiter=',', names=True, dtype=None, encoding='utf8')
    assert len(raw) == 4 * protocol['repetitions'] and len(summary) == 4
    names = ['Control_global_macro', 'Refined_global_macro_unstable', 'Refined_global_fine', 'Refined_local']
    rows = {str(r['case']): r for r in summary}
    assert set(rows) == set(names)
    for name in names:
        r = rows[name]
        trials = raw[raw['case'] == name]
        assert len(trials) == protocol['repetitions']
        assert np.array_equal(np.sort(trials['trial']), np.arange(protocol['repetitions']))
        for field, q in [('loop_median_ms', .5), ('loop_q25_ms', .25), ('loop_q75_ms', .75)]:
            assert np.isclose(r[field], np.quantile(trials['loop_ms'], q))
        assert np.isclose(r['setup_median_ms'], np.median(trials['setup_ms']))
        assert np.isclose(r['solver_total_median_ms'], np.median(trials['loop_ms'] + trials['setup_ms']))
        assert np.isfinite(trials['loop_ms']).all() and np.all(trials['loop_ms'] > 0)
        assert np.isclose(r['physical_duration_s'], r['dt_s'] * r['steps'])
        assert np.allclose(trials['physical_duration_s'], r['physical_duration_s'])
        if name == 'Refined_global_macro_unstable':
            assert r['steps'] == protocol['unstable_witness_steps']
            assert np.all(trials['amplitude_ratio'] > 1e5)
            assert np.isclose(r['physical_duration_s'], result['global_macro_divergence_time_s'])
        else:
            assert np.isclose(r['physical_duration_s'], result['physical_duration_s'])
            assert np.all(trials['amplitude_ratio'] < 10)
    for name, key in [('Refined_global_fine', 'global_fine_final_mass_l2_error'),
                      ('Refined_local', 'local_final_mass_l2_error')]:
        assert np.allclose(raw['relative_l2_error'][raw['case'] == name], result[key], rtol=1e-5)
    fine = rows['Refined_global_fine']['loop_median_ms']
    local = rows['Refined_local']['loop_median_ms']
    metrics = {'passed': True, 'speedup_local_vs_tested_global_fine': float(fine / local),
               'local_loop_time_reduction_vs_tested_global_fine': float(1-local/fine),
               'same_refined_mesh_and_physical_duration_for_speedup': True,
               'not_a_search_for_best_stable_global_time_step': True,
               'unstable_timing_excluded_from_speedup': True}
    (root / 'timing_analysis.json').write_text(json.dumps(metrics, indent=2)+'\n')
    labels = {'Control_global_macro': 'Témoin, global grand pas',
              'Refined_global_macro_unstable': 'Raffiné, global grand pas — instable',
              'Refined_global_fine': 'Raffiné, global petit pas', 'Refined_local': 'Raffiné, local'}
    lines = ['| Calcul | Durée physique | Préparation médiane | Propagation médiane | Intervalle interquartile |',
             '|---|---:|---:|---:|---:|']
    for name in names:
        r = rows[name]
        duration = f"{r['physical_duration_s']:.6f} s" if name.endswith('unstable') else '3,6 s'
        loop = f"{r['loop_median_ms']/1000:.4f} s" if not name.endswith('unstable') else f"{r['loop_median_ms']:.3f} ms"
        iqr = (f"{r['loop_q25_ms']/1000:.4f}–{r['loop_q75_ms']/1000:.4f} s" if not name.endswith('unstable') else
               f"{r['loop_q25_ms']:.3f}–{r['loop_q75_ms']:.3f} ms")
        lines.append(f"| {labels[name]} | {duration} | {r['setup_median_ms']:.3f} ms | {loop} | {iqr} |")
    (root/'timing_table.md').write_text('\n'.join(lines)+'\n')
    stable = [rows[n] for n in ['Control_global_macro', 'Refined_global_fine', 'Refined_local']]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), layout='constrained', width_ratios=[1.7, 1])
    median = np.array([r['loop_median_ms']/1000 for r in stable])
    low = median - np.array([r['loop_q25_ms']/1000 for r in stable])
    high = np.array([r['loop_q75_ms']/1000 for r in stable]) - median
    bars = axes[0].bar(['Témoin\nmaillage initial', 'Global petit pas\nmaillage raffiné', 'Local\nmaillage raffiné'], median,
                      color=['#b38b31', '#286b9c', '#078578'], yerr=[low, high], capsize=4)
    axes[0].bar_label(bars, labels=[f'{v:.3f} s' for v in median], padding=8)
    axes[0].set(ylabel='Boucle de propagation (s)', title='Calculs stables · 3,6 s physiques chacun', ylim=(0, max(median+high)*1.2))
    unstable = rows['Refined_global_macro_unstable']
    bars = axes[1].bar(['Global grand pas\nmaillage raffiné'], [unstable['loop_median_ms']], color='#b52e38',
                       yerr=[[unstable['loop_median_ms']-unstable['loop_q25_ms']],
                             [unstable['loop_q75_ms']-unstable['loop_median_ms']]], capsize=4)
    axes[1].bar_label(bars, labels=[f"{unstable['loop_median_ms']:.3f} ms"], padding=8)
    axes[1].set(ylabel='Boucle jusqu’au seuil (ms)',
                title=f"Témoin instable · {int(unstable['steps'])} pas seulement\nArrêt à {unstable['physical_duration_s']*1000:.3f} ms physiques",
                ylim=(0, unstable['loop_q75_ms']*1.25))
    for ax in axes:
        ax.grid(axis='y', alpha=.2)
    fig.suptitle(f"Temps CPU · {protocol['repetitions']} répétitions après {protocol['warmups']} échauffements", fontsize=16)
    fig.supxlabel('Kokkos Serial · un cœur CPU fixé · médiane et IQR · préparation, diagnostics et écritures exclus de la boucle', fontsize=10)
    fig.savefig(root/'04_temps_calcul.png', dpi=160)
    plt.close(fig)
    print('\n'.join(lines))
    print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', type=Path, default=Path(__file__).parent/'gmsh-inclusion-stability-results')
    analyze(p.parse_args().results)
