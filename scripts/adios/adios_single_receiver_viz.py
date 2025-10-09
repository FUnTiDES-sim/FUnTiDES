#!/usr/bin/env python
"""
Read and plot seismogram data from ADIOS BP5 file
Using the correct ADIOS2 Python API with FileReader
"""

import numpy as np
import matplotlib.pyplot as plt
from adios2 import FileReader

def read_adios_seismogram(filename):
    """
    Read seismogram data from ADIOS BP5 file

    Parameters:
    -----------
    filename : str
        Path to the BP5 file/directory

    Returns:
    --------
    data : numpy array
        Seismogram data
    """

    print(f"Opening file: {filename}")

    # Use FileReader to open the BP file
    with FileReader(filename) as bp_file:

        # List all available variables
        print("\nAvailable variables:")
        for var_name in bp_file.available_variables():
            print(f"  - {var_name}")

        # Try to read the AccousticReceiver variable
        var_name = "AccousticReceiver"

        # Inquire the variable to get its properties
        var = bp_file.inquire_variable(var_name)

        if var is not None:
            print(f"\nVariable '{var_name}' found:")
            print(f"  Shape: {var.shape()}")
            print(f"  Type: {var.type()}")

            # Read the entire variable
            data = bp_file.read(var_name)

            print(f"  Successfully read data with shape: {data.shape}")
            print(f"  Data range: [{np.min(data):.6f}, {np.max(data):.6f}]")

            return data
        else:
            print(f"\nVariable '{var_name}' not found in the file")
            print("Attempting to read with different variable names...")

            # Try alternative names
            alternative_names = ["AcousticReceiver", "acousticreceiver", "receiver"]

            for alt_name in alternative_names:
                var = bp_file.inquire_variable(alt_name)
                if var is not None:
                    print(f"  Found variable with name: '{alt_name}'")
                    data = bp_file.read(alt_name)
                    print(f"  Successfully read data with shape: {data.shape}")
                    return data

            # If still not found, try to read the first available variable
            available_vars = bp_file.available_variables()
            if available_vars:
                first_var = list(available_vars.keys())[0]
                print(f"\nReading first available variable: '{first_var}'")
                data = bp_file.read(first_var)
                print(f"  Successfully read data with shape: {data.shape}")
                return data

            print("No variables could be read from the file")
            return None

def plot_seismogram(data, title="Seismogram", receiver_index=0):
    """
    Plot the seismogram data

    Parameters:
    -----------
    data : numpy array
        Seismogram data (shape: n_receivers x n_timesteps)
    title : str
        Title for the plot
    receiver_index : int
        Which receiver to plot (if multiple receivers)
    """

    if data is None:
        print("No data to plot")
        return None

    # Handle different data shapes
    if data.ndim == 2:
        # If 2D, extract the specified receiver
        if data.shape[0] == 1:
            # Single receiver, multiple time steps
            trace = data[0, :]
            print(f"Plotting single receiver with {data.shape[1]} time samples")
        elif data.shape[1] == 1:
            # Multiple receivers, single time step (transpose)
            trace = data[:, 0]
            print(f"Plotting {data.shape[0]} samples")
        else:
            # Multiple receivers, multiple time steps
            trace = data[receiver_index, :]
            print(f"Plotting receiver {receiver_index} out of {data.shape[0]} receivers")
    else:
        # If 1D, plot directly
        trace = data.flatten()
        print(f"Plotting 1D data with {len(trace)} samples")

    # Create time axis
    n_samples = len(trace)
    time = np.arange(n_samples)

    # Create the plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

    # Time series plot
    ax1.plot(time, trace, 'b-', linewidth=0.5)
    ax1.set_xlabel('Time sample')
    ax1.set_ylabel('Amplitude')
    ax1.set_title(f'{title} - Time Series')
    ax1.grid(True, alpha=0.3)

    # Statistics
    print(f"\nSeismogram statistics:")
    print(f"  Number of samples: {n_samples}")
    print(f"  Min amplitude: {np.min(trace):.6f}")
    print(f"  Max amplitude: {np.max(trace):.6f}")
    print(f"  Mean amplitude: {np.mean(trace):.6f}")
    print(f"  Std deviation: {np.std(trace):.6f}")

    # Amplitude spectrum
    try:
        from scipy import signal
        frequencies, psd = signal.periodogram(trace)
        ax2.semilogy(frequencies[:len(frequencies)//2], psd[:len(psd)//2])
        ax2.set_xlabel('Normalized Frequency')
        ax2.set_ylabel('Power Spectral Density')
        ax2.set_title('Frequency Content')
        ax2.grid(True, alpha=0.3)
    except ImportError:
        ax2.text(0.5, 0.5, 'scipy not available for frequency analysis',
                ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('Frequency Content (scipy required)')

    plt.tight_layout()
    plt.show()

    return trace

def save_seismogram_to_text(data, output_filename="seismogram.txt"):
    """
    Save seismogram data to a text file

    Parameters:
    -----------
    data : numpy array
        Seismogram data
    output_filename : str
        Output filename
    """
    if data is None:
        print("No data to save")
        return

    # Flatten the data for saving
    flat_data = data.flatten()

    # Save with metadata
    header = f"Seismogram data\n"
    header += f"Original shape: {data.shape}\n"
    header += f"Total samples: {len(flat_data)}\n"
    header += f"Min: {np.min(flat_data):.6e}, Max: {np.max(flat_data):.6e}"

    np.savetxt(output_filename, flat_data, fmt='%.6e', header=header)
    print(f"\nSeismogram saved to '{output_filename}'")

# Main execution
if __name__ == "__main__":
    import sys

    # Get file path from command line or use default
    if len(sys.argv) > 1:
        bp_file = sys.argv[1]
    else:
        bp_file = "receivers.bp"

    print(f"Reading ADIOS file: {bp_file}")
    print("=" * 50)

    # Read the seismogram data
    sismo_data = read_adios_seismogram(bp_file)

    if sismo_data is not None:
        # Plot the seismogram
        plot_seismogram(sismo_data, title="Acoustic Receiver Seismogram")

        # Optional: Save to text file
        save_seismogram_to_text(sismo_data, "seismogram_output.txt")

        # Optional: Save to numpy format for later use
        np.save("seismogram_data.npy", sismo_data)
        print(f"Seismogram data saved to 'seismogram_data.npy'")
    else:
        print("\nFailed to read seismogram data from file")
