Use "pip install -r requirements.txt" to install requirements package to your virtual environment (Py == 3.8).

DFROStNET:
- Python version must be the same as the one used to train the AI model. Since I used 3.8 to train the DFROStNET file that currently exists (as of April 14th 2025), we must use 3.8 to generate predictions.
python == 3.8

PtyChoPy compile files:
- the C files must be compiled where you're running the UI (if you don't, you'll get an import error for the call_ wrappers).
- Compiled files can only be run on the same Py --v it was installed on. So you may have to recompile:
    1. Ensure you have MS Visual Studio for C++ installed and pathed in your system environments (make sure to select C++ build tools during installation).
    2. Navigate to PtyChoPy/Cython_wrappers folder in cmd.
        2.1. Navigate into the wrapper folder (e.g., fc_decalage_temporelle_vecteur_t_de_t0_BEN)
        2.2. Run "python setup.py build_ext --inplace"
        2.3. It will compile and output a modified '.pyd'
    3. Repeat the above for each _BEN folder in Cython_wrappers.
    4. For fc_projection_experimentale:
        4.1. _BEN is out of date. Installing the FFTW in C is difficult to reproduce but was necessary in MATLAB.
        4.2. Instead, we use _SYD which implements the pyFFTW implementation to perform the same computation in roughly the equivalent time.
        4.3. If you can manage to setup _BEN properly, that's probably still a bit faster but I can't figure it out so good luck. 
- if you move / rename files, YOU MUST RECOMPILE. This must be done for each of the _BEN files in Cython_wrappers. 


KNOWN BUGS:
if you get:
- an import error on the _BEN files in the PtyChoPy: rebuild the compiled files
- after creating a custom cmap, the first time you switch back to a standard map for some reason
    the cmap generator still appears. if you close it, it carries on as normal. 
    No clue why, but it's not stopping the rest from working so it's not a priority.
