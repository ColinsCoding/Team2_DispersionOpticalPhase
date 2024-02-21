README

Data is generated using VPI (Virtual Photonics Inc.) code for binary WDM four wavelength communication. 

Data format is intensity modulation (IM). Fiber propagation length is 35km of SMF fiber. 


In the subfolder dispersion_only, the fiber only has dispersion. 
In the subfolder dispersion_and_nonlinearity, the fiber has both dispersion and nonlinearity. 
 
Data is saved in 8 csv files. RX_lambda1.csv is the received data in wavelength 1 and TX_lambda1.csv is the transmit (labels) in wavelength 1.

Data has only one row, which is the samples of the received signal. Data is 25 Gbit/s for each wavelength channel, sampling rate is 400 GS/s. There are 16384 symbols in each channel in total.

For Coherent (QAM16) data, we are interested in the short propagation fiber distance (100m) and the main non-idealities come from the locol oscillator linewidth and the phase imbalance from the Hybrid module in coherent detection. 

16QAM_25G_ideal contains data where the linewidth is 0 and the imbalance is 0.
16QAM_25G_IQ_imbalance_*_deg contains data where the Hybrid module has IQ imbalance.
16QAM_25G_LO_linewidth_*MHz contains data where the local oscillator laser has linewidth. 

For QAM16 data, the saved TX are transmitted symbols (00,01,10,11 for IQ channel respectively), the saved RX are the decoded signal at 25GHz.
