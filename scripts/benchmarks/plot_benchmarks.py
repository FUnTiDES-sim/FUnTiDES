import os
import re
import json
import argparse
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple


def parse_args():
    p = argparse.ArgumentParser(description="Plot Python vs C++ benchmark results.")
    p.add_argument('--python-dir',
                   required=True,
                   help='Directory containing Python benchmark JSON files (tX.json).')
    p.add_argument('--cpp-dir',
                   required=True,
                   help='Directory containing C++ benchmark JSON files (tX.json).')
    p.add_argument('--save',
                   default='./benchmark_comparison.png',
                   help='Path to save the resulting plot (PNG).')
    return p.parse_args()


def extract_data(folder: str, extractor) -> Tuple[List[int], Dict[str, Dict[int, float]]]:
    """
    Scan folder for JSON files with thread count pattern and extract benchmark data.
    This function searches for JSON files matching the pattern `t<threads>.json` 
    (e.g., t1.json, my_test_t16.json) and extracts benchmark performance data 
    using the provided extractor function.
    Parameters
    ----------
    folder : str
        Path to the directory containing JSON benchmark files.
    extractor : callable
        Function to extract benchmark data from JSON files. Should accept two 
        parameters: file path (str) and a dictionary to store extracted data.
    Returns
    -------
    per_thread_data : Dict[int, Dict[str, Dict[int, float]]]
        Nested dictionary containing benchmark data organized by thread count.
        Structure: {threads: {benchmark_name: {order: time, ...}, ...}, ...}
        Expected benchmark names:
        - 'structured_FE_init'
        - 'structured_FE_one_step'
        - 'unstructured_FE_init'
        - 'unstructured_FE_one_step'
    Notes
    -----
    - Files that don't match the pattern `t<threads>.json` are skipped
    - Parse errors for individual files are caught and printed as warnings
    """
    thread_pattern = re.compile(r'.*t(\d+)\.json$')
    per_thread_data: Dict[int, Dict[str, Dict[int, float]]] = {}

    for entry in os.scandir(folder):
        if entry.is_file() and entry.name.endswith('.json'):
            m = thread_pattern.match(entry.name)
            if not m:
                continue
            threads = int(m.group(1))

            if threads not in per_thread_data:
                per_thread_data[threads] = {
                    'structured_FE_init': {},
                    'structured_FE_one_step': {},
                    'unstructured_FE_init': {},
                    'unstructured_FE_one_step': {}
                }

            try:
                print(f"Processing {entry.path} for {threads} threads...")
                extractor(entry.path, per_thread_data[threads])
            except Exception as e:
                print(f"Warning: failed to parse {entry.path}: {e}")

    return per_thread_data


def extract_python_benchmarks(json_file, results):

    with open(json_file, 'r') as f:
        data = json.load(f)

        for benchmark in data.get('benchmarks', []):
            name = benchmark.get('name', '')
            params = benchmark.get('params', {})

            # Extract order and struct parameters (assume size of list > 0)
            if 'struct' in params and isinstance(params['struct'], list):
                order = params['struct'][0]
                is_structured = True
                on_nodes = params['struct'][-1]
            elif 'unstruct' in params and isinstance(params['unstruct'], list):
                order = params['unstruct'][0]
                is_structured = False
                on_nodes = params['unstruct'][-1]
            else:
                raise ValueError("Missing 'struct' or 'unstruct' parameter for mapping")

            # Determine benchmark type
            is_init = 'solver_fe_init' in name
            is_one_step = 'solver_one_step' in name

            # TODO for now skip onNodes
            if on_nodes:
                continue

            # Map to appropriate category (in ms)
            time_ms = benchmark['stats']['mean'] * 1000
            fill_data(results,
                      is_structured,
                      is_init, is_one_step,
                      order,
                      time_ms)


def extract_cpp_benchmarks(json_file, results):

    with open(json_file, 'r') as f:
        data = json.load(f)

        for benchmark in data.get('benchmarks', []):
            name = benchmark.get('name', '')
            label = benchmark.get('label', '')

            # Determine if structured or unstructured
            is_structured = 'SolverStructFixture' in name

            # Determine benchmark type
            is_init = 'FEInit' in name
            is_one_step = 'OneStep' in name

            # Extract order from label
            if label:
                match = re.search(r'Order=(\d+)', label)
                order = int(match.group(1))
                on_nodes = 'OnNodes=1' in label
            else:
                raise ValueError("Missing 'label' parameter for mapping")

            # TODO for now skip onNodes
            if on_nodes:
                continue

            # Map to appropriate category (in ms)
            time_ms = benchmark['real_time']
            fill_data(results,
                      is_structured,
                      is_init, is_one_step,
                      order,
                      time_ms)


def fill_data(data: Dict[str, Dict[int, float]],
              is_structured: bool,
              is_init: bool,
              is_one_step: bool,
              order: int,
              time_ms: float):

    # Check exclusivity of init/one_step flags
    if is_init and is_one_step:
        raise ValueError("Cannot be both init and one_step")
    if not is_init and not is_one_step:
        raise ValueError("Must be either init or one_step")

    # Check order is strictly positive
    if order <= 0:
        raise ValueError("Order must be a positive integer")

    if is_structured and is_init:
        data['structured_FE_init'][order] = time_ms
    elif is_structured and is_one_step:
        data['structured_FE_one_step'][order] = time_ms
    elif not is_structured and is_init:
        data['unstructured_FE_init'][order] = time_ms
    elif not is_structured and is_one_step:
        data['unstructured_FE_one_step'][order] = time_ms


def plot_benchmarks(python_data: Dict[int, Dict[str, Dict[int, float]]],
                    cpp_data: Dict[int, Dict[str, Dict[int, float]]],
                    save_path: str = None):
    """
    Plot benchmark comparison results for Python and C++ implementations.

    Creates a multi-panel figure with one subplot per benchmark category, showing
    runtime performance across different thread counts and polynomial orders.

    Parameters
    ----------
    python_data : Dict[int, Dict[str, Dict[int, float]]]
        Nested dictionary containing Python benchmark results.
        Structure: {threads: {category: {order: runtime_ms}}}
    cpp_data : Dict[int, Dict[str, Dict[int, float]]]
        Nested dictionary containing C++ benchmark results.
        Structure: {threads: {category: {order: runtime_ms}}}
    save_path : str, optional
        File path to save the figure. If None, displays the plot interactively.
        Default is None.

    Notes
    -----
    The function generates a figure with 4 subplots, one for each category:
    - Structured FE Init
    - Structured FE One Step
    - Unstructured FE Init
    - Unstructured FE One Step

    Each subplot displays:
    - X-axis: Thread counts
    - Y-axis: Runtime in milliseconds
    - Multiple curves representing different polynomial orders
    - Python implementations shown with solid lines and circle markers
    - C++ implementations shown with dashed lines and x markers
    - Different colors distinguish polynomial orders (using tab10/tab20 colormaps)

    The figure includes a legend for each subplot and a main title.
    If save_path is provided, saves the figure at 180 DPI.
    """
    categories = [
        ('structured_FE_init', 'Structured FE Init'),
        ('structured_FE_one_step', 'Structured FE One Step'),
        ('unstructured_FE_init', 'Unstructured FE Init'),
        ('unstructured_FE_one_step', 'Unstructured FE One Step')
    ]

    # All thread counts (sorted)
    thread_counts = sorted(set(python_data.keys()) | set(cpp_data.keys()))

    # One subplot per category
    fig, axs = plt.subplots(1, len(categories), figsize=(4 * len(categories), 4), squeeze=False)
    axs = axs[0]

    for ax, (cat_key, cat_title) in zip(axs, categories):
        # Collect all orders present for this category
        orders = set()
        for t in thread_counts:
            orders.update(python_data.get(t, {}).get(cat_key, {}).keys())
            orders.update(cpp_data.get(t, {}).get(cat_key, {}).keys())
        orders = sorted(o for o in orders if o is not None)

        if not orders:
            ax.set_title(f"{cat_title}\n(no data)")
            ax.axis('off')
            continue

        cmap_name = 'tab10' if len(orders) <= 10 else 'tab20'
        cmap = plt.get_cmap(cmap_name, len(orders))

        for idx, order in enumerate(orders):
            color = cmap(idx)

            # Python series
            py_vals = []
            py_threads = []
            for t in thread_counts:
                v = python_data.get(t, {}).get(cat_key, {}).get(order)
                if v is not None:
                    py_threads.append(t)
                    py_vals.append(v)

            # C++ series
            cpp_vals = []
            cpp_threads = []
            for t in thread_counts:
                v = cpp_data.get(t, {}).get(cat_key, {}).get(order)
                if v is not None:
                    cpp_threads.append(t)
                    cpp_vals.append(v)

            if py_vals:
                ax.plot(py_threads, py_vals,
                        marker='o', markersize=4, linestyle='-',
                        color=color, markerfacecolor='none',
                        label=f"Py O{order}")
            if cpp_vals:
                ax.plot(cpp_threads, cpp_vals,
                        marker='x', linestyle='--',
                        color=color,
                        label=f"C++ O{order}")

        ax.set_title(cat_title)
        ax.set_xlabel("Threads")
        ax.set_ylabel("Runtime (ms)")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xticks(thread_counts)
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        ax.grid(alpha=0.3, which='both')
        ax.legend(fontsize='small', ncol=2)

    fig.suptitle('Benchmark Comparison by Category (Orders as curves)', fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    if save_path:
        plt.savefig(save_path, dpi=180)
        print(f"Saved figure to {save_path}")
    else:
        plt.show()
    plt.close(fig)


def main():
    args = parse_args()

    python_data = extract_data(args.python_dir, extract_python_benchmarks)
    cpp_data = extract_data(args.cpp_dir, extract_cpp_benchmarks)

    if len(python_data) != len(cpp_data):
        print("Warning: Python and C++ sets differ.")

    plot_benchmarks(python_data, cpp_data, save_path=args.save)


if __name__ == "__main__":
    main()
