function [status,cmdout] = proteowizard(inputfilename,outputfilepath)
% PROTEOWIZARD: This function converts a raw file to a mzML or other XML
%     file using Proteowizard's msconvert function
%
% Syntax:
% [status,cmdout] = proteowizard(inputfilename,outputfilepath)
%
% Input:
% inputfilename: string. Full path of input expt. data file.
% outputfilepath: string. Folder where produced file is saved.
%
% Output:
% status: nonzero integer. Command exit status. 0 means successful
% cmdout: string. Output of operating system command
%
% Note:
% It is important to maintain Proteowizard to its newest version.
% When Proteowizard is updated, we recommend uninstall the old version
%     ASAP. The path to Proteowizard is semi-hard coded, older version has
%     priority unless this file is changed or the older version has been
%     uninstalled.
%
% Examples:
% status = proteowizard('a filename','a filepath')
%
% See also:
% GETMSDATA  GETINDEXEDMZMLDATA  GETMZMLDATA  GETMZXMLDATA
%

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 08/18/2023

currentLoc = pwd;
listing = dir();
proteowizardpath = '';
if any(cellfun(@(x) contains (x,'proteowizardloc.txt'),{listing.name}))
    fileID=fopen('proteowizardloc.txt')
    newLoc = fscanf(fileID,'%c');
    fclose(fileID);
    listing=dir(newLoc);
    if any(cellfun(@(x) contains (x,'msconvert.exe'),{listing.name}))
        proteowizardpath=newLoc;
    end
    if strcmpi(proteowizardpath(end),'\') || strcmpi(proteowizardpath(end),'/')
        proteowizardpath(end) = '';
    end
else
    newLoc = [currentLoc,'\ProgramFiles\ProteoWizard 3.0.22353.3f760db'];
    listing=dir(newLoc);
    if any(cellfun(@(x) contains (x,'msconvert.exe'),{listing.name}))
        proteowizardpath=newLoc;
    end
end

if ~exist(proteowizardpath,'dir')
    [~,proteowizardpath,~] = uigetfile('*.exe','Please locate folder with Proteowizard "msconvert" executable');
    fileID = fopen('proteowizardloc.txt','w');
    formatSpec = '%c';
    fprintf(fileID,formatSpec,proteowizardpath);
    fclose(fileID);
end

% Default folder path to proteowizard. If folder does not exist, program
%     will ask user to specify a new one
%{
%This program overwrites the getproteowizardloc file to your personal path
if ~exist(proteowizardpath,'dir')
    [~,proteowizardpath,~] = uigetfile('*.exe','Please Locate "msconvert" Program');
    if ispc
        currentloc_split = strsplit(mfilename('fullpath'),'\');
        proteowizardlocfile = strjoin([currentloc_split(1:end-1),'getproteowizardloc.m'],'\');
    else
        currentloc_split = strsplit(mfilename('fullpath'),'/');
        proteowizardlocfile = strjoin([currentloc_split(1:end-1),'getproteowizardloc.m'],'/');
    end
    fileID = fopen(proteowizardlocfile,'w');
    newfilecontent = 'function proteowizardpath = getproteowizardloc()\nproteowizardpath = ''%s'';\nend';
    fprintf(fileID,newfilecontent,proteowizardpath);
    fclose(fileID);
%}


% Navigate to Proteowizard installation folder to run DOS command
cd(proteowizardpath)
if strcmpi(outputfilepath(end),'\') || strcmpi(outputfilepath(end),'/')
    outputfilepath(end) = '';
end

% DOS command
file_type = 'mzML'; % Type of file to convert to; default is mzML
message = strcat('Converting .RAW file to .',file_type);
g = msgbox(message);
peakPickingOutput = 'true [2,3]'; % PeakPicking; default is true [2,3]
command = strcat('msconvert',{' "'},inputfilename,{'" -o "'},outputfilepath,{'" --'},file_type,{' --filter "peakPicking '},peakPickingOutput,{'"'});
[status,cmdout] = dos(command{1});
close(g);

cd(currentLoc);
end
