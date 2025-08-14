import numpy as np
import tensorflow as tf
import scipy
import os 
import matplotlib.pyplot as plt
import threading

parameters = {'Default_DFROStNET': [[1, 3, 5], 3, 0.21, np.sqrt(1 / (10**(7 / 10))), 5, 575],
              'dummy_load_tester': [[7, 5, 3, 1], 3, 0.25, 0.1, 1, 512],
              'best_train':        [[3, 5, 5, 5], 1, 0.18, 0.35, 4, 166],
              'best_val':          [[5, 5, 7, 3], 1, 0.53, 0.5, 3, 1860]}


"""Ptychography thread to make fasters. """
class PtychoThread(threading.Thread):
    """Thread used to reduce lag for ptycho constraint."""
    def __init__(self, data, callback=None):
        self.trace, self.pulse, self.time, self.delay = data
        self.callback = callback  # Optional callback to notify when processing is complete

    def run(self):
        self.calc_switch = self.apply_ptychographic_constraint(self.time, self.delay, self.trace, self.pulse)
        if self.callback:
            self.callback(self.calc_switch)

    """ Functions for ptychographic constraint. 
        Grandfathered in from Ben's (Benoit Brizzard, 2019, I think) MATLAB C-wrapped functions. 
        I took the steps to perform the constraint and reproduced them in Py for this application 
        (Note that Py is much slower, this is done for simplicity because the constraint is applied only once). """
    def apply_ptychographic_constraint(self, time, delay, trace, pulse, true_switch = None):
        og_shape = len(delay)
        new_shape = [1051, 1501]

        delay = self.interp_vector(delay, new_shape[0])
        time = self.interp_vector(time, new_shape[1])
        trace = self.interpolate_batch(trace, new_shape[1], new_shape[0])
        pulse = tf.complex(self.interpolate_batch(tf.math.real(pulse), 0, new_shape[1]), self.interpolate_batch(tf.math.imag(pulse), 0, new_shape[1]))
        N, N0 = len(time), len(delay)
        dt, dt0 = time[1] - time[0], delay[1] - delay[0]

        # # step 1 get e field matrix from trace
        trace_t = self.batchwise_2D_ifft_(trace)

        # step 2 create shifted e field matrix
        M_decalage_t0 = self.fc_decalage_temporelle_batch_tf(N0, N, dt0, dt, pulse)

        # step 3 somme sur decales
        M_decalage_t0_sum = self.fc_somme_sur_les_decales_batch_tf(N0, N, dt0, dt, M_decalage_t0)

        # # step 4 equation
        pred_switch = self.sum_colonne_normale_batch(M_decalage_t0_sum, trace_t, N0, N)

        # interpolate back down to # of points
        pred_switch = self.interpolate_batch(pred_switch, 0, og_shape)

        # uses label to determine phase shift if simulation data
        if true_switch is not None:
            pred_switch = self.roll_and_pad(pred_switch, true_switch)

        pred_switch = self.normalize_switch(np.abs(pred_switch))

        return pred_switch
    
    @staticmethod
    def interp_vector(vector, new_len):
        new_shape = (new_len / vector.shape[0]) if vector.ndim == 1 else (new_len[0] / vector.shape[0], new_len[1] / vector.shape[1])
        return scipy.ndimage.zoom(vector, new_shape)
        
    @staticmethod
    def batchwise_2D_ifft_(trace):
        dummy_trace = tf.cast(trace[:, :, :, 0], tf.complex64)
        dummy_trace = tf.transpose(tf.signal.fftshift(dummy_trace, axes = 1), perm = [0, 2, 1])  # Swap -2 (axis 2) with the last axis
        trace_t = tf.signal.ifft(dummy_trace)
        trace_t = tf.signal.ifftshift(tf.transpose(trace_t, perm = [0, 2, 1]), axes = 1)
        trace_t = tf.expand_dims(trace_t, axis = -1)
        return trace_t
    
    @staticmethod
    def roll_and_pad(labels, recovereds):
        signal_length = tf.shape(labels)[1]
        batch_size = tf.shape(labels)[0]

        # Convert inputs to complex64 for FFT operations
        labels = tf.cast(labels, tf.complex64)
        recovereds = tf.cast(recovereds, tf.complex64)

        # Compute FFT of labels and recovereds
        label_fft = tf.signal.fftshift(tf.signal.fft(tf.signal.ifftshift(labels, axes=1)), axes=1)
        recovered_fft = tf.signal.fftshift(tf.signal.fft(tf.signal.ifftshift(recovereds, axes=1)), axes=1)

        # Cross-spectrum calculation
        cross_spectrum = label_fft * tf.math.conj(recovered_fft)

        # Compute the phase shift
        phase_shift = tf.math.angle(tf.signal.ifft(tf.signal.fftshift(cross_spectrum, axes=1)))

        # Find the shift for each signal in the batch
        shift_indices = tf.cast(tf.argmax(phase_shift, axis = 1), dtype = tf.int32)  # Shape: (batch_size,)
        shifts = -(shift_indices - signal_length // 2 ) // 152 # division is required to scale down shifts (trial and error, larger division gives smaller error)
        
        # Compute the shifted indices for each signal
        range_indices = tf.range(signal_length, dtype = tf.int32)
        shifted_indices = range_indices[None, :] - shifts[:, None]
        clipped_indices = tf.clip_by_value(shifted_indices, 0, signal_length - 1)
        aligned_signals = tf.gather(recovereds, clipped_indices, batch_dims=1)
        return aligned_signals

    @staticmethod
    def fc_decalage_temporelle_batch_tf(N0, N, dt0, dt, V_t_batch):
        """
        Optimized TensorFlow version for handling a batch of vectors.
        
        Parameters:
            N0 (int): Number of rows in the output matrices.
            N (int): Number of columns in the output matrices.
            dt0 (float): Time step for the rows.
            dt (float): Time step for the columns.
            V_t_batch (Tensor): Tensor of input vectors of shape (batch_size, vector_length).
        
        Returns:
            Tensor: Tensor of output matrices of shape (batch_size, N0, N).
        """
        batch_size, vector_length = tf.shape(V_t_batch)[0], tf.shape(V_t_batch)[1]        
        N_dt0_s_dt = dt0 / dt
        N_pts_plus_decalage = tf.cast((N0 - 1) * N_dt0_s_dt, tf.int32)
        moitie = N_pts_plus_decalage // 2

        # Create row indices for each vector in the batch
        row_indices = tf.range(N0, dtype=tf.float32) * dt0 / dt
        row_indices = tf.cast(tf.round(row_indices), tf.int32)

        # Prepare indices for the columns
        col_indices = tf.range(-moitie, N - moitie, dtype=tf.int32)

        # Create a 2D grid of indices for rows and columns
        col_grid = tf.expand_dims(col_indices, axis=0) + tf.expand_dims(row_indices, axis=1)

        # Apply boundary conditions
        col_grid_clipped = tf.clip_by_value(col_grid, 0, vector_length - 1)

        # Gather values for each batch
        def process_vector(v):
            return tf.gather(v, col_grid_clipped, axis=0)

        # Process each vector in the batch
        M_decalage_t0_batch = tf.map_fn(process_vector, V_t_batch, fn_output_signature = tf.complex64)

        M_decalage_t0_batch = tf.image.flip_left_right(M_decalage_t0_batch)

        return M_decalage_t0_batch

    @staticmethod
    def fc_somme_sur_les_decales_batch_tf(N0, N, dt0, dt, M_decalage_t0):
        batch_size = tf.shape(M_decalage_t0)[0]
        N_dt0_s_dt = dt0 / dt
        N_pts_plus_somme = int(tf.floor((N_dt0_s_dt + 1.0) / 2.0))
        extended_N = N + 2 * N_pts_plus_somme - (0 if N_dt0_s_dt % 2 else 1)

        # Initialize the extended tensor
        M_t_t0_loc = tf.zeros((batch_size, N0, extended_N), dtype=M_decalage_t0.dtype)

        # Fill the initial columns
        initial_columns = tf.repeat(
            tf.expand_dims(M_decalage_t0[:, :, 0], axis=-1), 
            repeats=N_pts_plus_somme - (0 if N_dt0_s_dt % 2 else 1), 
            axis=-1
        )
        M_t_t0_loc = tf.concat([initial_columns, M_t_t0_loc[:, :, tf.shape(initial_columns)[-1]:]], axis=-1)

        # Fill the main columns
        main_columns = M_decalage_t0
        M_t_t0_loc = tf.concat(
            [
                M_t_t0_loc[:, :, :N_pts_plus_somme - 1],
                main_columns,
                M_t_t0_loc[:, :, N_pts_plus_somme - 1 + tf.shape(main_columns)[-1]:],
            ],
            axis=-1,
        )

        # Fill the last columns
        last_columns = tf.repeat(
            tf.expand_dims(M_decalage_t0[:, :, -1], axis=-1), 
            repeats=N_pts_plus_somme, 
            axis=-1
        )
        M_t_t0_loc = tf.concat(
            [
                M_t_t0_loc[:, :, :N_pts_plus_somme - 1 + tf.shape(main_columns)[-1]],
                last_columns,
            ],
            axis=-1,
        )

        # Compute the sum for the shifted matrices
        M_decalage_t0_sum = tf.zeros((batch_size, N0, N), dtype=M_decalage_t0.dtype)
        for i in range(int(N_dt0_s_dt)):
            shifted = M_t_t0_loc[:, :, i : i + N]
            M_decalage_t0_sum += shifted

        return M_decalage_t0_sum

    @staticmethod
    def interpolate_batch(batch_tensor, target_height, target_width):
        if batch_tensor.ndim > 3:
            # Use tf.image.resize for resizing the batch
            resized_tensor = tf.image.resize(
                batch_tensor, 
                size=(target_height, target_width), 
                method="bicubic"  # or "nearest", "bicubic", etc., depending on the desired interpolation method
            )
        else:
            expanded_tensor = tf.expand_dims(batch_tensor, axis=-1)  # (batch_size, length, 1)
            expanded_tensor = tf.expand_dims(expanded_tensor, axis=1)  # (batch_size, 1, length, 1)
            resized_tensor = tf.image.resize(
                expanded_tensor, 
                size=(1, target_width),  # Resize only along the length axis
                method="bicubic"
            )

            # Squeeze back to remove extra dimensions
            resized_tensor = tf.squeeze(resized_tensor, axis=[1, -1])  # (batch_size, target_length)
        return resized_tensor

    @staticmethod
    def sum_colonne_normale_batch(M_Obj, M_trace_t, N0, N):
        M_Obj_summed = tf.reduce_sum(M_Obj, axis=1)  # Shape: (batch_size, N)
        M_Obj_summed_conj = tf.math.conj(M_Obj_summed)  # Shape: (batch_size, N)
        M_Obj_summed_conj = tf.cast(M_Obj_summed_conj, tf.complex64)
        numerator = tf.einsum('bn,bnij->bni', M_Obj_summed_conj, M_trace_t)  # Shape: (batch_size, N0, N)
        numerator = tf.transpose(numerator, perm=[0, 2, 1])
        denominator = tf.reduce_sum(tf.abs(M_Obj_summed) ** 2, axis=1)  # Shape: (batch_size,)
        Obj_t = tf.abs(tf.reduce_sum(numerator, axis = 2) / tf.expand_dims(tf.cast(denominator, tf.complex64), axis = 1))  # Shape: (batch_size, N0)
        max_values = tf.reduce_max(tf.abs(Obj_t), axis = 1, keepdims = True)  # Shape: (batch_size, 1)
        Obj_t = Obj_t / max_values  # Normalize
        return Obj_t

    @staticmethod
    def normalize_switch(tensor):
        tensor_min = tf.reduce_min(tensor, axis = -1, keepdims = True)
        tensor_max = tf.reduce_max(tensor, axis = -1, keepdims = True)
        normalized_tensor = (tensor - tensor_min) / (tensor_max - tensor_min)
        return normalized_tensor
    
    @staticmethod
    def fc_decalage_temporelle_Matrice_t_de_t0_SYD(N0, N, dt0, dt, V_t):
        N_dt0_s_dt = dt0 / dt
        N_pts_plus_decalage = int((N0 - 1) * N_dt0_s_dt)
        moitie = N_pts_plus_decalage // 2

        M_decalage_t0 = np.zeros((N0, N), dtype = V_t.dtype)
        for row in range(N0):
            ind = int(row * dt0 / dt) + 1
            for col in range(moitie, moitie + N):
                # First range: Before the shift index
                if col < ind - 1:
                    M_decalage_t0[row, col - moitie] = V_t[0]
                # Second range: Between the shift indices
                elif ind - 1 <= col < ind + N - 2:
                    M_decalage_t0[row, col - moitie] = V_t[col - ind + 1]
                # # Third range: After the shift index
                else:
                    M_decalage_t0[row, col - moitie] = V_t[N - 1]
        return M_decalage_t0

    @staticmethod
    def fc_somme_sur_les_dt0_decales_Matrice_de_t_decalee_de_t0_BEN(N0, N, dt0, dt, M_decalage_t0):
        N_dt0_s_dt = dt0 / dt
        N_pts_plus_somme = int(np.floor((N_dt0_s_dt + 1.0) / 2.0)) #int((np.floor((N_dt0_s_dt + 1.0) / 2.0)))
        M_decalage_t0_sum = np.zeros((N0, N), dtype = np.complex128)
        if N_dt0_s_dt % 2:  # Odd case
            M_t_t0_loc = np.zeros((N0, N + 2 * N_pts_plus_somme - 1), dtype=np.complex128)

            # Fill the initial columns
            for j in range(N_pts_plus_somme - 1):
                M_t_t0_loc[:, j] = M_decalage_t0[:, 0]

            # Fill the main columns
            for j in range(N):
                M_t_t0_loc[:, N_pts_plus_somme - 1 + j] = M_decalage_t0[:, j]

            # Fill the last columns
            for j in range(N_pts_plus_somme):
                M_t_t0_loc[:, N_pts_plus_somme - 1 + N + j] = M_decalage_t0[:, -1]

        else:  # Even case
            M_t_t0_loc = np.zeros((N0, N + 2 * N_pts_plus_somme), dtype=np.complex128)
            # Fill the initial columns
            for j in range(N_pts_plus_somme):
                M_t_t0_loc[:, j] = M_decalage_t0[:, 0]

            # Fill the main columns
            for j in range(N):
                M_t_t0_loc[:, N_pts_plus_somme + j] = M_decalage_t0[:, j]

            # Fill the last columns
            for j in range(N_pts_plus_somme):
                M_t_t0_loc[:, N_pts_plus_somme + N + j] = M_decalage_t0[:, -1]

        # Add shifted matrices
        for indice_de_N_dt0_s_dt in range(int(N_dt0_s_dt)):
            for j in range(N):
                M_decalage_t0_sum[:, j] += M_t_t0_loc[:, indice_de_N_dt0_s_dt + j]

        return M_decalage_t0_sum

    @staticmethod
    def sum_colonne_normale(M_Obj, M_trace_t, N0, N):   
        M_Obj = np.sum(M_Obj, axis = 0)
        numerator = np.zeros((N0, N), dtype = np.complex128)
        for i in range(N0):
            numerator[i, :] = np.conj(M_Obj) * M_trace_t[:, i]
        denominator = np.sum(np.abs(M_Obj) ** 2)
        Obj_t = np.sum(numerator, axis = 1) / denominator
        Obj_t /= np.max(Obj_t)
        return Obj_t

""" Class for main DFROStNET functionalities. It will be called from Master window as necessary. """
class DFROStNET(threading.Thread):
    def __init__(self, model_path = None, callback = None):
        super().__init__()
        # load necessary files from ai training
        # self.model_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Models", "Default_DFROStNET.weights.h5") if model_path is None else model_path
        # self.model_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Models", "dummy_load_tester.weights.h5") if model_path is None else model_path
        self.model_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Models", "best_train.weights.h5") if model_path is None else model_path
        # self.model_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Models", "best_val.weights.h5") if model_path is None else model_path

        self.points = 128  # Corresponds to number of neurons in the architecture
        self.callback = callback
        self.model = None
        self.model_name = None
        self.load_model(self.model_path)

    def get_model_parameters(self):
        filename = os.path.basename(self.model_path)
        model_name = filename.split('.')[0]
        self.model_params = parameters[model_name]

    def run(self, data, estimate_spectrum):
        self.raw_trace, self.raw_delay, self.raw_wvl, self.raw_freq, self.raw_time = map(data.get, ['trace', 'delay', 'wvl', 'freq', 'time'])
        self.est_spectrum = estimate_spectrum
        self.make_prediction()
    
    def make_prediction(self):
        # proceed with prediction, as a threaded calculation
        self.generate_pulse()
        self.generate_switch()

    """ UI functions: for init, plotting, etc. """
    def load_model(self, path):
        self.model_path = path
        try:
            # replace call parameters if using a different model
            self.get_model_parameters()
            self.model = self.create_MR_model('multires', *self.model_params)
            self.model.load_weights(self.model_path)

        except Exception as e:
            print(f"Error loading model: {e}")

    def create_MR_block(self, input_matrix, filters, MR_kernels, conv_kernel):
        input_shape = input_matrix.shape

        # MR convolutions
        filtered = []
        for i in range(len(MR_kernels)):
            stride = (input_shape[1] - MR_kernels[i] + (2*MR_kernels[i]//2)) // (input_shape[1] - 1)
            x = tf.keras.layers.Conv2D(filters = filters, 
                                    kernel_size = (MR_kernels[i], MR_kernels[i]), 
                                    strides = (stride, stride), 
                                    padding = 'same', activation = 'relu')(input_matrix)
            filtered.append(x)
        output_matrix = tf.keras.layers.Concatenate(axis = -1)(filtered)

        # post MR convolution
        stride = (output_matrix.shape[1] - conv_kernel + 2*(conv_kernel//2)) // (output_matrix.shape[1]//2 - 1)
        output_matrix = tf.keras.layers.Conv2D(filters = output_matrix.shape[-1] * 2, 
                                            kernel_size = (conv_kernel, conv_kernel), 
                                            strides = (stride, stride), 
                                            padding = 'same', activation = 'relu')(output_matrix)
        return output_matrix

    def create_MR_model(self, name, MR_kernels = [5, 3, 1], conv_kernel = 1, dropout = 0, stddev = 0, dense_layers = 1, dense_neurons = 512):
        input_shape = (self.points, self.points, 1)
        inputs = tf.keras.Input(shape = input_shape)
        wgn_layer = tf.keras.layers.GaussianNoise(stddev)(inputs)

        # multi res blocks
        multires_output1 = self.create_MR_block(wgn_layer, input_shape[-1], MR_kernels, conv_kernel)
        multires_output2 = self.create_MR_block(multires_output1, multires_output1.shape[-1], MR_kernels, conv_kernel)
        multires_output3 = self.create_MR_block(multires_output2, multires_output2.shape[-1], MR_kernels, conv_kernel)
        
        # flatten, drop
        x = tf.keras.layers.Flatten()(multires_output3)
        x = tf.keras.layers.Dropout(dropout)(x)

        # dense and output layers
        for _ in range(dense_layers):
            x = tf.keras.layers.Dense(dense_neurons, activation = 'relu')(x)

        switch_output = tf.keras.layers.Dense(3*self.points, activation = 'linear')(x)
        phi_output = tf.keras.layers.Dense(self.points, activation = 'linear')(x)
        model = tf.keras.Model(inputs = inputs, outputs = [phi_output, switch_output], name = name)
        return model

    """ Functional standalone methods used throughout. """
    @staticmethod
    def interpolate_data(to_interp, length = 128):
        new_shape = (length[0] / to_interp.shape[0], length[1] / to_interp.shape[1]) if to_interp.ndim == 2 else (length / to_interp.shape[0],)
        return scipy.ndimage.zoom(to_interp, new_shape, order = 3)
        
    @staticmethod
    def _ifft_(spectrum, axis = 1):
        return tf.signal.ifftshift(tf.signal.ifft(tf.signal.fftshift(spectrum, axes = axis)), axes = axis)
    
    @staticmethod
    def _fft_(pulse, axis = 1):
        return tf.signal.fftshift(tf.signal.fft(tf.signal.ifftshift(pulse, axes = axis)), axes = axis)
    
    """ Required for AI predictions. """
    def interpolate_raw_dataset(self, points = [128, 128]):
        interp_trace = self.interpolate_data(self.raw_trace, points)
        interp_trace /= np.max(np.abs(interp_trace))
        interp_delay = self.interpolate_data(self.raw_delay, points[1])
        interp_wvl = self.interpolate_data(self.raw_wvl, points[0])
        interp_time = self.interpolate_data(self.raw_time, points[0])
        return interp_trace, interp_delay, interp_wvl, interp_time
    
    @staticmethod
    def normalize_switch(tensor):
        tensor_min = tf.reduce_min(tensor, axis = 1, keepdims = True)
        tensor_max = tf.reduce_max(tensor, axis = 1, keepdims = True)
        normalized_tensor = (tensor - tensor_min) / (tensor_max - tensor_min)
        return normalized_tensor
    
    def generate_pulse(self):
        # interpolate down to 128, 128 points first
        self.interp_trace, interp_delay, self.interp_wvl, self.pred_time = self.interpolate_raw_dataset()
        self.interp_trace = tf.reshape(self.interp_trace, (1, *self.interp_trace.shape, 1))
        if self.est_spectrum.shape[-1] != self.interp_trace.shape[1]:   
            self.est_spectrum = self.interpolate_data(self.est_spectrum, self.interp_trace.shape[1])

        # generate ai prediction: self.pred_pulse, self.pred_spectrum, self.approx_switch, self.approx_trace
        pred_phase, pred_switch = self.model.predict(self.interp_trace, verbose = 0)
        self.pred_pulse = self.return_prediction(self.est_spectrum, pred_phase*np.pi)
        pred_switch = self.interpolate_data(pred_switch[0], 4*self.points).reshape(1, -1)
        self.approx_switch = np.flip(self.normalize_switch(pred_switch), axis = 1)
        self.approx_trace = self.FROStNET(self.pred_pulse, self.approx_switch)

        # send to plot in main_win
        self.callback({
            "pulse": self.pred_pulse[0] if self.pred_pulse.ndim == 2 else self.pred_pulse,
            "time": self.pred_time,
            "wvl": self.interp_wvl,
            "pro_trace": self.interp_trace[0, :, :, 0].numpy() if self.interp_trace.ndim == 4 else self.interp_trace.numpy(),
            "trace": self.approx_trace[0, :, :, 0].numpy() if self.approx_trace.ndim == 4 else self.approx_trace.numpy(), 
            "switch": self.approx_switch[0].numpy() if self.approx_switch.ndim == 2 else self.approx_switch.numpy(),
            "delay": interp_delay,
        }, "DFROStNET")

    def return_prediction(self, amp, phase):
        spectrum = self.recombine_complex_vector(amp, phase) 
        spectrum /= np.max(spectrum, axis = -1, keepdims = True)
        pulse = self._ifft_(spectrum)
        pulse -= np.abs(np.mean(np.concatenate([pulse[:, :10], pulse[:, -10:]])))
        pulse *= np.hanning(pulse.shape[-1])
        pulse /= np.max(pulse, axis = -1, keepdims = True)
        return pulse
    
    @staticmethod
    def recombine_complex_vector(amplitude, phase):
        amplitude_complex = tf.cast(amplitude, tf.complex64)
        phase_complex = tf.complex(tf.zeros_like(phase), phase)
        return amplitude_complex * tf.exp(phase_complex)

    def generate_switch(self):
        # generate switch prediction: self.pred_switch
        self.ptycho_thread = PtychoThread([self.interp_trace, self.pred_pulse, self.raw_time, self.raw_delay], 
                            callback = self.handle_switch_signal)
        self.ptycho_thread.run()

    def handle_switch_signal(self, switch):
        switch = scipy.ndimage.gaussian_filter1d(np.abs(switch), sigma = 10) 
        switch = switch - (np.max(switch[:]) - np.mean(switch[:, :10])) - np.mean(switch[:, -10:])
        self.pred_switch = self.normalize_switch(self.interpolate_data(switch[0], 4*self.points).reshape(1, -1))
        self.pred_trace = self.FROStNET(self.pred_pulse, self.pred_switch)
        self.pred_switch = self.interpolate_data(self.pred_switch[0], self.points)
        self.pred_delay = self.interpolate_data(self.raw_delay, self.points)

        self.callback({
            "pulse": self.pred_pulse[0] if self.pred_pulse.ndim == 2 else self.pred_pulse,
            "time": self.pred_time,
            "wvl": self.interp_wvl,
            "pro_trace": self.interp_trace[0, :, :, 0].numpy() if self.interp_trace.ndim == 4 else self.interp_trace.numpy(),
            "trace": self.pred_trace[0, :, :, 0].numpy() if self.pred_trace.ndim == 4 else self.pred_trace.numpy(), 
            "switch": self.pred_switch[0] if self.pred_switch.ndim == 2 else self.pred_switch,
            "delay": self.pred_delay,
        }, "DFROStNET")

    def FROStNET(self, pulse, switch):
        points = pulse.shape[-1]
        switch_shifted = np.lib.stride_tricks.sliding_window_view(switch, window_shape = (points, ) , axis = 1)
        product = switch_shifted * pulse[:, np.newaxis, :]    
        trace = np.abs(self._fft_(product, axis = -1))**2
        trace = trace.reshape(trace.shape[0], trace.shape[1], trace.shape[2], -1)
        trace /= np.max(trace, axis = (1, 2), keepdims = True)
        trace = tf.transpose(trace, perm=[0, 2, 1, 3])
        trace = tf.image.resize(trace, (points, points))
        return trace
    
if __name__ == "__main__":

    dfn = DFROStNET()
