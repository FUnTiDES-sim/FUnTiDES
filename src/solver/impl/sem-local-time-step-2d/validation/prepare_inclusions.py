"""Rasterize named inclusions and preserve independent time-step targets.

The target steps are design inputs, not instructions to the current 2:1 solver.
No rounding to an integer macro-step ratio is performed here.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def positive(value):
    assert isinstance(value, (int, float)) and math.isfinite(value) and value > 0, value
    return value


def contains(inc, x, z):
    cx, cz = inc["center"]
    if inc["shape"] == "box":
        sx, sz = inc["size"]
        return abs(x - cx) < sx / 2 and abs(z - cz) < sz / 2
    if inc["shape"] == "disk":
        return (x - cx) ** 2 + (z - cz) ** 2 < inc["radius"] ** 2
    raise ValueError(f"Unsupported shape: {inc['shape']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    here = Path(__file__).parent
    parser.add_argument("--config", type=Path, default=here / "inclusion_models.json")
    parser.add_argument("--output", type=Path, default=here / "inclusion-results")
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    mesh = config["mesh"]
    # The experimental executable reuses the existing Q2 validation mesh.
    assert mesh == {"lx": 6000, "lz": 3000, "ex": 60, "ez": 30, "order": 2}
    ex, ez, lx, lz = (mesh[k] for k in ("ex", "ez", "lx", "lz"))
    pulse = config["initial_field"]
    assert 0 <= pulse["x"] <= lx and 0 <= pulse["z"] <= lz
    positive(pulse["sigma"])
    positive(config["duration_s"])
    positive(config["partition"]["macro_dt_s"])
    positive(config["partition"]["cfl"])
    names = [case["name"] for case in config["cases"]]
    assert len(set(names)) == len(names)
    for case in config["cases"]:
        name = case["name"]
        assert name and all(c.isalnum() or c in "_-" for c in name)
        regions = [config["background"], *case["inclusions"]]
        for region in regions:
            for key in ("vp", "rho", "target_dt_s"):
                positive(region[key])
        for inc in case["inclusions"]:
            cx, cz = inc["center"]
            assert math.isfinite(cx) and math.isfinite(cz)
            if inc["shape"] == "box":
                sx, sz = map(positive, inc["size"])
            elif inc["shape"] == "disk":
                sx = sz = 2 * positive(inc["radius"])
            else:
                raise ValueError(inc["shape"])
            assert 0 <= cx - sx / 2 < cx + sx / 2 <= lx
            assert 0 <= cz - sz / 2 < cz + sz / 2 <= lz
        counts = [0] * len(regions)
        rows = []
        for k in range(ez):
            for i in range(ex):
                x, z = (i + .5) * lx / ex, (k + .5) * lz / ez
                matches = [j + 1 for j, inc in enumerate(case["inclusions"]) if contains(inc, x, z)]
                if len(matches) > 1:
                    raise ValueError(f"Overlapping inclusions in element {i + k * ex}")
                region_id = matches[0] if matches else 0
                region = regions[region_id]
                counts[region_id] += 1
                rows.append((i + k * ex, region["vp"], region["rho"], region_id, region["target_dt_s"]))
        assert all(counts), f"Empty region in {name}: {counts}"
        output = args.output / name
        output.mkdir(parents=True, exist_ok=True)
        with (output / "model.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["element", "vp", "rho", "region_id", "target_dt_s"])
            writer.writerows(rows)
        p = config["partition"]
        boundary = case.get("boundary_conditions", config.get("boundary_conditions", "rigid"))
        modes = {"rigid": 0, "absorbing": 1, "free_top_absorbing_sides": 2, "free_top_rigid_sides": 3}
        assert boundary in modes
        free_top = boundary in ("free_top_absorbing_sides", "free_top_rigid_sides")
        receiver_depth = config.get("receivers", {}).get("depth_m", 0)
        assert math.isfinite(receiver_depth) and 0 <= receiver_depth <= lz
        assert abs(receiver_depth / 50 - round(receiver_depth / 50)) < 1e-8
        if free_top:
            assert pulse["z"] < lz, "The pressure excitation must be below the pressure-free surface"
            assert receiver_depth > 0, "Pressure receivers exactly on the free surface would measure zero"
        reference_precision = config.get("reference_time_accumulation", "float32")
        assert reference_precision in ("float32", "float64")
        (output / "parameters.txt").write_text(
            f"{pulse['x']} {pulse['z']} {pulse['sigma']} {config['duration_s']}\n"
            f"{p['macro_dt_s']} {p['cfl']}\n{modes[boundary]}\n{int(reference_precision == 'float64')}\n{receiver_depth}\n")
        metadata = {"name": name, "mesh": mesh, "initial_field": pulse,
                    "boundary_conditions": boundary,
                    "receivers": {"depth_m": receiver_depth, "z": lz-receiver_depth,
                                  "x_min": 500, "x_max": 5500, "count": 51},
                    "initial_condition_profile": "gaussian_dirichlet_image" if free_top else "gaussian",
                    "reference_time_accumulation": reference_precision,
                    "duration_s": config["duration_s"], "partition": p,
                    "regions": [{"id": j, "elements": counts[j], **r} for j, r in enumerate(regions)],
                    "target_steps_are_executed": False,
                    "actual_integrator": "LF-LTS 2:1; fixed projection from velocity diagnostics"}
        (output / "model.json").write_text(json.dumps(metadata, indent=2) + "\n")
        print(f"{name}: region element counts {counts}; targets {[r['target_dt_s'] for r in regions]} s")
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "models.json").write_text(json.dumps(config, indent=2) + "\n")


if __name__ == "__main__":
    main()
