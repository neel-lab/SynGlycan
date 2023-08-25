function installglycopat
% INSTALLGLYCOPAT function is to install the GLYCOPAT. It consists of
% several pre-check subroutines and installation steps.
%
% 1. System check and display basic information
% 2. Add the GLYCOPAT directory to the Path
% 3. Add the java libraries to the java class path
% 4. Check if GLYCOPAT works properly

% Author: Gang Liu
% Date Lastly Updated: 12/01/14

clc;clear;  %clear command window and workspace

fprintf(1,'Installation starts \n');
fprintf(1,'-------------------------------------------------------------\n');
fprintf(1,'1. Display Matlab system information.\n');
dispSystemInfo();
checkVersions();

fprintf(1,'\n2. GlycoPAT Installation\n');
fprintf(1,'\t2.1 Set up environment path\n');
addGlycoPATPath();

%fprintf(1,'\t2.2 Check Proteowizard folder path\n');
%fprintf(1,'\t\tPlease specify Proteowizard folder path\n');
%getproteowizard();
%fprintf(1,'\t\tProteowizard folder is set\n');

fprintf(1,'\n3.Installation check\n');
fprintf(1,'\t 3.1.Folder path check\n');
functiontest();
fprintf(1,'---------------------------------------------------------------\n');
fprintf(1,'Installation ends.\n');
synglycangui
end

function functiontest()
%functiontest check path setup and main folder existence
dir1_GLYCOPATToolbox = 'ProgramFiles';
isDirsGLYCOPATToolboxExist = exist(dir1_GLYCOPATToolbox,'dir');
if (isDirsGLYCOPATToolboxExist)
    fprintf(1,'\t\tValidation:GlycoPAT is installed properly\n');
else
    error('GlycoPAT:Installation:INCORRECTPATH','\t\tGLYCOPAT is not installed properly\n');
end
end

function addGlycoPATPath()
modulenames = {'ProgramFiles'};
rootloc = which('installSynGlycan.m');
[f,~,~] = fileparts(rootloc);
addpath(f);
for i = 1:length(modulenames)
    glycopatpath = fullfile(f,modulenames{i});
    addpath(genpath(glycopatpath));
    s = savepath;
    if (s ~= 0)
        warning('GlycoPAT:Installation:PathNotAdded',...
            '\nGLYCOPAT directorie and subdirectories were not added to the Pathdef.m. ');
    else
        fprintf(1,'\t\t Installation successful\n');
    end
end
end

function dispSystemInfo()
%systemDisplay display matlab information
%
%  See also COMPUTER,ISUNIX,IS32BIT
curMatlab = ver('matlab');
currMatlabVer = curMatlab.Version;

fprintf(1,'  \tMatlab platform: %s version %s\n',getSystemsByte,currMatlabVer);

if(ispc)
    sysComputer='  PC\n';
elseif(ismac)
    sysComputer=   ' Mac\n';
else
    sysComputer=  '  Unix/Linux\n';
end

fprintf(1,sprintf('\tPC system: %s',sysComputer));
end

function result = is32bit()
%IS32BIT True for the 32bit  version of MATLAB.
%   IS32BIT returns 1 for 32bit  versions of MATLAB and 0 otherwise.
%
%   See also COMPUTER, ISUNIX, ISMAC,ISPC.

result = strcmp(computer,'PCWIN') || strcmp(computer,'GLNX86');
end

function systemByte = getSystemsByte()
%getSystemsByte return a string for the MATLAB bytes
%   getSystemsByte returns 32bit and 64bit based on current Matlab platform.
%
%   See also is32bit, ISUNIX, ISMAcomputC,ISPC.
if(is32bit())
    systemByte = '32 bit';
else
    systemByte = '64 bit';
end
end

function isSuggestedVersion=checkVersions( )
% checkVersions True for the recommended Matlab version or higher
%    checkVersions returns 1 for suggested Matlab verion and 0 otherwise
%
%  See also verLessThan
StableVersionForGlycoToolbox = '8.2';

curMatlab = ver('matlab');
currMatlabVer = curMatlab.Version;
% toolboxParts = getParts(toolboxver(1).Version);
% verParts = getParts(verstr);
isSuggestedVersion=~verLessThan('matlab',StableVersionForGlycoToolbox);
if(~isSuggestedVersion)
    warning('GLYCOPAT:Installation:MATLABVERSIONOWER',sprintf('Matlab Version %s or higher is recommended',StableVersionForGlycoToolbox));
else
    sprintf('Current version %s is stable for GPAT',currMatlabVer);
end

end

function getproteowizard()
% get ProteoWizard path and save to ""
[~,proteowizardpath,~] = uigetfile('*.exe','Please Locate "msconvert" Program');
if ispc
    currentloc_split = strsplit(mfilename('fullpath'),'\');
    proteowizardlocfile = strjoin([currentloc_split(1:end-1),'Function\xmlreader\getproteowizardloc.m'],'\');
else
    currentloc_split = strsplit(mfilename('fullpath'),'/');
    proteowizardlocfile = strjoin([currentloc_split(1:end-1),'Function/xmlreader/getproteowizardloc.m'],'/');
end
fileID = fopen(proteowizardlocfile,'w');
newfilecontent = 'function proteowizardpath = getproteowizardloc()\nproteowizardpath = ''%s'';\nend';
fprintf(fileID,newfilecontent,proteowizardpath);
fclose(fileID);
end
