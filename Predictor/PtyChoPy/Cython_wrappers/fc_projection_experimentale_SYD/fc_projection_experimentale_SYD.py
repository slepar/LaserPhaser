import numpy as np
import pyfftw
import os 

def ifftshift1d(input_array):
    """Equivalent to the C function ifftshift1D."""
    n = input_array.shape[0]
    pivot = (n // 2) if n % 2 == 0 else (n - 1) // 2
    return np.concatenate((input_array[pivot:], input_array[:pivot]))

def fftshift1d(input_array):
    """Equivalent to the C function fftshift1D."""
    n = input_array.shape[0]
    pivot = (n // 2) if n % 2 == 0 else (n + 1) // 2
    return np.concatenate((input_array[pivot:], input_array[:pivot]))

def fc_projection_experimentale_py(m_produit_t_t0_it, trace_xp_int):
    """
    Optimized Python equivalent of the C function fc_projection_experimentale_BEN.
    Inputs:
        m_produit_t_t0_it: Complex input matrix.
        trace_xp_int: Experimental intensity data.
    Outputs:
        Updated complex matrix, experimental amplitude, and error.
    """
    # Validate inputs
    if m_produit_t_t0_it.shape != trace_xp_int.shape:
        raise ValueError("Input dimensions do not match.")
    
    M, N = m_produit_t_t0_it.shape[0], m_produit_t_t0_it.shape[1]

    # Prepare FFTW arrays and plans
    fft_in = pyfftw.empty_aligned(N, dtype="complex128")
    fft_out = pyfftw.empty_aligned(N, dtype="complex128")
    fft = pyfftw.FFTW(fft_in, fft_out, direction="FFTW_FORWARD")
    ifft = pyfftw.FFTW(fft_out, fft_in, direction="FFTW_BACKWARD")

    # Initialize outputs
    m_produit_t_t0_it_plus_1 = np.zeros_like(m_produit_t_t0_it, dtype="complex128")
    trace_amp_omg_t0_it = np.zeros_like(m_produit_t_t0_it, dtype="complex128")

    # Scalars for error calculation
    t_xp_carre, produit_xp_recons, t_recons_carre = 0.0, 0.0, 0.0

    for row in range(M):
        # Load the row of data into the input array
        slice = m_produit_t_t0_it[row, :] / np.max(np.abs(m_produit_t_t0_it[row, :]))
        fft_in[:] = ifftshift1d(slice)
        fft()
        shifted_fft = fftshift1d(fft_out)

        # Calculate reconstruction intensity and apply experimental constraints
        trace_recons_int = np.abs(shifted_fft)**2
        t_xp_carre += np.sum(trace_xp_int[row, :]**2)
        produit_xp_recons += np.sum(trace_xp_int[row, :] * trace_recons_int)
        t_recons_carre += np.sum(trace_recons_int**2)

        # Update amplitudes based on constraints
        scaling_factors = np.sqrt(trace_xp_int[row, :] / trace_recons_int) 
        scaling_factors = np.nan_to_num(scaling_factors)  # Replace NaN/Inf with valid values
        shifted_fft *= scaling_factors

        # Store amplitude values
        trace_amp_omg_t0_it[row, :] = shifted_fft

        # Perform IFFT and normalization
        fft_out[:] = ifftshift1d(shifted_fft)
        ifft()
        m_produit_t_t0_it_plus_1[row, :] = fftshift1d(fft_in) / N

    # Normalize updated matrix
    div_num = np.mean(np.abs(m_produit_t_t0_it_plus_1)) 
    m_produit_t_t0_it_plus_1 /= div_num

    # Compute error
    err = np.sqrt(1 - (produit_xp_recons**2) / (t_xp_carre * t_recons_carre))
    return m_produit_t_t0_it_plus_1, trace_amp_omg_t0_it, err

# Example usage
if __name__ == "__main__":

    path = r"C:\Users\lepar\OneDrive - INRS\Documents\INRS\MSc\thesis\PtyPy\recovered_traces\20201201_FRost_3mm fused silica.txt"

    recovered_data = {}

    for filename in os.listdir(path):
        if filename.endswith(".txt"):
            name = os.path.splitext(filename)[0]  # Get the filename without the .txt extension
            filepath = os.path.join(path, filename)
            try:
                recovered_data[name] = np.loadtxt(filepath)
            except ValueError:      # handles the case where pulse is complex
                recovered_data[name] = np.loadtxt(filepath, dtype = complex)

    # print(recovered_data.keys())

    # print(recovered_data['processed_trace'].shape, recovered_data['processed_freq'].shape, recovered_data['processed_delay'].shape )

    trace_int = (recovered_data['processed_trace'].T)**2
    time_trace = np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(trace_int)))
    # plt.subplot(211)
    # plt.imshow(np.abs(trace))
    # plt.subplot(212)
    # plt.imshow(np.abs(time_trace))
    # plt.show()

    updated_matrix, experimental_amplitude, error = fc_projection_experimentale_py(time_trace, trace_int)

    # plt.subplot(221)
    # plt.imshow(np.abs(trace))
    # plt.subplot(222)
    # plt.imshow(np.abs(time_trace))
    # plt.subplot(223)
    # plt.imshow(np.abs(updated_matrix))
    # plt.subplot(224)
    # plt.imshow(np.abs(experimental_amplitude))

    # plt.show()
