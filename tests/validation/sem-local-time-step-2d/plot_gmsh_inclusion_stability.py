"""Verify and display the four-run inclusion stability test."""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/funtides-lts-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Results and reports of the validations
RESULTS = Path(__file__).resolve().parents[3] / "documentation" / "acoustic" / "local-time-step"


def plot(root):
    result = json.loads((root / 'test_result.json').read_text())
    calibration = json.loads((root / 'calibration.json').read_text())
    checks = json.loads((root / 'mesh_checks.json').read_text())
    assert checks['outside_mesh_identical_to_baseline_at_1e_minus6_m']
    assert result['passed'] and result['fine_outside_inclusion'] == 0
    assert result['first_growth_in_inclusion']
    assert 2800 <= result['first_growth_x_m'] <= 3200 and 700 <= result['first_growth_z_up_m'] <= 1100
    assert result['global_macro_amplification_at_stop'] > 1e5
    assert max(result['control_max_amplification'], result['local_max_amplification'],
               result['global_fine_max_amplification']) < 10
    assert result['max_local_energy_balance_error'] < 1e-4
    assert result['reference_delta_l2'] < .001
    assert max(result['local_final_mass_l2_error'], result['global_fine_final_mass_l2_error'],
               result['local_receiver_relative_l2_error']) < .01
    data = np.genfromtxt(root / 'stability_history.csv', delimiter=',', names=True)
    traces = np.genfromtxt(root / 'receivers.csv', delimiter=',', names=True)
    assert len(data) == result['macro_steps'] + 1
    assert np.isclose(data['time'][-1], result['physical_duration_s'])
    for name in data.dtype.names:
        if name != 'global_macro_peak':
            assert np.isfinite(data[name]).all(), name
    witness = np.isfinite(data['global_macro_peak'])
    assert np.array_equal(witness, np.arange(len(data)) <= np.flatnonzero(witness)[-1])
    assert np.isclose(data['time'][witness][-1], result['global_macro_divergence_time_s'])
    assert np.all(np.diff(data['dissipation']) >= -1e-12)
    e0 = data['energy'][0]
    balance = np.abs(data['energy'] + data['dissipation'] - e0) / e0
    assert np.allclose(balance, data['balance_error'], rtol=1e-4, atol=1e-10)
    assert np.isclose(balance.max(), result['max_local_energy_balance_error'], atol=1e-10)
    shape = (len(data), 51)
    reference = traces['reference'].reshape(shape)
    local = traces['local'].reshape(shape)
    fine = traces['global_fine'].reshape(shape)
    assert np.isfinite(reference).all() and np.isfinite(local).all() and np.isfinite(fine).all()
    assert np.allclose(traces['time'].reshape(shape), data['time'][:, None])
    err = np.linalg.norm(local - reference) / np.linalg.norm(reference)
    assert np.isclose(err, result['local_receiver_relative_l2_error'])
    dt = result['macro_dt_s'] * 1000
    small = result['fine_dt_s'] * 1000
    colors = {'global_macro_peak': '#b52e38', 'global_fine_peak': '#286b9c',
              'local_peak': '#078578', 'control_peak': '#b38b31'}
    labels = {'global_macro_peak': f'Global raffiné, grand pas {dt:.4f} ms',
              'global_fine_peak': f'Global raffiné, petit pas {small:.4f} ms',
              'local_peak': f'Local, {dt:.4f} / {small:.4f} ms',
              'control_peak': f'Témoin sans sur-raffinement, {dt:.4f} ms'}
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), layout='constrained')
    for ax in axes[0]:
        for name in colors:
            ax.semilogy(data['time'], data[name] / data[name][0], label=labels[name], color=colors[name], lw=1.4)
        ax.axhline(1e5, color='#b52e38', ls=':', lw=1)
        ax.set(xlabel='Temps physique (s)', ylabel='Amplitude maximale / amplitude initiale')
        ax.grid(alpha=.2)
    axes[0, 0].set_title('Le grand pas global diverge ; les trois témoins restent bornés')
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].set(xlim=(0, result['global_macro_divergence_time_s'] * 1.1), title='Zoom sur le déclenchement dans l’inclusion')
    axes[0, 1].annotate(f"Premier facteur 10 à {result['first_tenfold_growth_time_s']:.5f} s\n"
                        f"x = {result['first_growth_x_m']:.2f} m, profondeur = {3000-result['first_growth_z_up_m']:.2f} m",
                        xy=(result['first_tenfold_growth_time_s'], 10), xytext=(.04, .85), textcoords='axes fraction',
                        fontsize=9, arrowprops=dict(arrowstyle='->'),
                        bbox=dict(facecolor='white', edgecolor='none', alpha=.95))
    axes[1, 0].plot(data['time'], data['energy'] / e0, label='Énergie locale / E₀')
    axes[1, 0].plot(data['time'], data['dissipation'] / e0, label='Dissipation cumulée / E₀')
    axes[1, 0].set(xlabel='Temps physique (s)', title='Bilan local avec absorption')
    axes[1, 0].legend(fontsize=9)
    axes[1, 1].semilogy(data['time'], np.maximum(data['balance_error'], 1e-16), label='Erreur du bilan énergétique')
    axes[1, 1].semilogy(data['time'], np.maximum(data['local_l2_error'], 1e-16), label='Erreur locale du champ vs référence')
    axes[1, 1].set(xlabel='Temps physique (s)', ylabel='Erreur relative', title='Précision et dissipation contrôlées')
    axes[1, 1].legend(fontsize=9)
    fig.suptitle('Petit pas nécessaire dans l’inclusion à 600 m/s : test réussi', fontsize=17)
    fig.supxlabel('Le témoin divergent est arrêté au facteur 10⁵ ; les trois calculs stables continuent jusqu’à 3,6 s.', fontsize=10)
    fig.savefig(root / '01_necessite_petit_pas.png', dpi=160)
    plt.close(fig)

    gain = np.minimum(data['time'] ** 2, 20)[:, None]
    extent = (.5, 5.5, data['time'][-1], 0)
    amplitude = np.max(np.abs(reference * gain))
    fig, axes = plt.subplots(1, 3, figsize=(14, 6), layout='constrained')
    for ax, values, label in zip(axes, (reference, local, local-reference), ('Référence', 'Local', 'Local − référence')):
        im = ax.imshow(values * gain, extent=extent, aspect='auto', cmap='RdBu_r',
                       vmin=-amplitude, vmax=amplitude, interpolation='nearest')
        ax.set(xlabel='x du récepteur (km)', ylabel='Temps (s)', title=label)
    fig.colorbar(im, ax=axes, label='Pression × gain · même échelle')
    fig.suptitle(f'Sismogrammes à 50 m de profondeur · erreur L² sans gain : {err:.4%}', fontsize=15)
    fig.supxlabel('Gain commun min[(t / 1 s)², 20] · aucune normalisation par trace', fontsize=10)
    fig.savefig(root / '02_sismogrammes.png', dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 7), layout='constrained')
    xs = traces['x'].reshape(shape)[0]
    for ax, x in zip(axes.flat, [1000, 2200, 3000, 4000]):
        j = np.argmin(abs(xs-x))
        ax.plot(data['time'], reference[:, j] * gain[:, 0], label='Référence', color='#286b9c')
        ax.plot(data['time'], local[:, j] * gain[:, 0], label='Local', color='#078578', ls='--')
        ax.plot(data['time'], (local-reference)[:, j] * gain[:, 0], label='Erreur', color='#555555', lw=.8)
        ax.set(title=f'x = {xs[j]/1000:g} km', xlabel='Temps (s)', ylabel='Pression × gain')
        ax.grid(alpha=.2)
    axes[0, 0].legend()
    fig.suptitle('Traces locales et référence : même gain et mêmes axes')
    fig.savefig(root / '03_traces.png', dpi=160)
    plt.close(fig)
    (root / 'analysis_checks.json').write_text(json.dumps({'passed': True, 'raw_csv_consistency': True,
        'outside_mesh_identical': True, 'receiver_relative_l2_error_recomputed': float(err),
        'spectral_estimates_are_not_certified_bounds': True}, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', type=Path, default=RESULTS / 'gmsh-inclusion-stability-results')
    plot(p.parse_args().results)
