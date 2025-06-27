# SubSeaGasLeak

SubSeaGasLeak is a python model able to simulate the rising of a plume of gas bubble trough a water column.

## Installation

To install the model, clone the repository, and install the requirements in `requirements.txt` using the command `pip install -r requirements.txt`.

## Run the model

### With a graphical interface

To run the model, use the command `python run_bubble.py`. In the graphical interface, the simulation results folder can be opened with the button `Open results folder`.

### Without a graphical interface

To run the model without a graphical interface, use the command `python run_bubble.py --no_gui --sim_ID SIM_ID`, with a file named `SIM_ID_raw_request.json` in the `requests` folder (It will be created if it doesn't exist when the model is run with the GUI at least once, otherwise it must be create manually). `SIM_ID` is a number between 0 and 99999, that must be 5 digits (eg. 00123 and not 123). The results will be in the `results/SIM_ID` folder.

## Acknowledgements

The development of this software has been co-funded by the European Commission, DG-ECHO in the framework of the MANIFESTS Genius project (grant agreement n° 101140390).
