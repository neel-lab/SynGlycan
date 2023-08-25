function [proteowizardpath] = getproteowizardloc()
currentLoc = pwd;
newLoc = [currentLoc,'\ProgramFiles\ProteoWizard 3.0.22353.3f760db'];
listing=dir(newLoc);

if any(cellfun(@(x) contains (x,'msconvert.exe'),{listing.name}))
    proteowizardpath=newLoc;
else
    [~,proteowizardpath,~] = uigetfile('*.exe','Please locate folder with Proteowizard "msconvert" executable');
end

%{
lastbsloc=strfind(currentLoc,'\');
lastbsloc=lastbsloc(end);
pathCell = strsplit(currentLoc,'\');
onPath = strcmpi("SynGlycan", pathCell);
isonPath = any(onPath);
cellNum=find(onPath,1);
pathCell=pathCell(1:cellNum);
synGlycanPath=strjoin(pathCell,'\');
if extractAfter(currentLoc,lastbsloc) == "SynGlycan"
    proteowizardpath = [currentLoc,'\ProteoWizard 3.0.22353.3f760db'];
elseif isonPath
    proteowizardpath = [synGlycanPath,'\ProteoWizard 3.0.22353.3f760db'];
else
    [~,proteowizardpath,~] = uigetfile('*.exe','Please locate folder with Proteowizard "msconvert" executable');
    pathCell = strsplit(proteowizardpath,'\');
    onPath = strcmpi("SynGlycan", pathCell);
    cellNum=find(onPath,1);
    pathCell=pathCell(1:cellNum);
    synGlycanPath=strjoin(pathCell,'\');
proteowizardpath = strrep(currentLoc,'\Function\xmlreader','\ProteoWizard 3.0.22353.3f760db')
C:\Users\neel\Box\cbe-neel\cbe-neel-shared\Software\SynGlycan\ProgramFiles\ProteoWizard 3.0.22353.3f760db
 during compilation
 data_file_path = fullfile(ctfroot, 'ProteoWizard 3.0.22353.3f760db');
%}

end