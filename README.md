# SynGlycan
## MS glycoproteomics software customized for SynGlycan toolkit
This the code for GlycoPAT2 that is customized for the analysis of SynGlycan probes.\
To use the program, you can download the code and install it at your desired location.
Subsequently change directory to that location in the MATLAB environment.\
Write "install SynGlycan". All that does is add the downloaded files to the MATLAB path
Then double click the SynGlycan.mlapp to run the program.\
Please note that it is important for you to start the program from the installation folder, otherwise the relative dependencies may be missed. Additionally, you will need to download proteowizard independently, in order to convert .RAW files to .mat files.\
The program can take .RAW and .mzML files as input. These are converted into a _QA.mat file (essentially a dataframe) that is used for subsequent analysis.
