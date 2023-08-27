# SynGlycan
## MS glycoproteomics software customized for SynGlycan toolkit
1. This the code for GlycoPAT2 that is customized for the analysis of SynGlycan probes.<br>
2. To use the program, you can download the code to your desired location (say c:\User\software\).
Subsequently change directory to that installation directory in the MATLAB environment. The program is best used with MATLAB 2022b or later, but it could work with some previous versions also <br>
3. Write "installSynGlycan". This runs the associated .m file. All this does is add the downloaded files to the MATLAB path
Then double click the synglycangui.mlapp to run the program.<br>
4. Please note that it is important for you to start the program from the installation folder, otherwise the relative dependencies may be missed. Additionally, you will need to download proteowizard independently, in order to convert .RAW files to .mat files.<br>
5. The program can take .RAW and .mzML files as input. These are converted into a _QA.mat file (essentially a dataframe) that is used for subsequent analysis. Once you make the _QA.mat file once, it is better to use that as a starting material for all subsequent analysis rather than starting with .RAW or .mzmL and repeating the process all over again.<br>
Good luck and email us in case you need more help!
