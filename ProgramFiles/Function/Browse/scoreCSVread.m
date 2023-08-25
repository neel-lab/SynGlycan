function [headerpar,msdatatype,datatyp,scoredata] = scoreCSVread(fname)
% SCORECSVREAD: Read data from scoring result CSV file in RFC 4180 format
%
% Syntax:
%   [headerpar,msdatatype,datatyp,scoredata] = scoreCSVread(fname)
%
% Input:
% fname: .csv result file full name(including file path)
%
% Output:
% headerpar: structure, .csv file header. 
%     Currently fields are:
%         pepfilename: digestion file name.
%         DataDirectory: (legacy) when experiment data is in DTA format,
%             this is the folder DTA files were stored.
%         xmlfilename: full name of .mat or .mzXML file, including path.
%         fragMod: fragmentation modes that were analyzed.
%         MS1tolUnit: MS1 tolerence unit.
%         MS1tol: MS1 tolerence value.
%         MS2tolUnit: MS2 tolerence unit for each fragmentation mode.
%         MS2tol: MS2 tolerence value for each fragmentation mode.
%         npFrag: max. number of peptide bond cleavage allowed for each
%             fragmentation mode.
%         ngFrag: max. number of glycan bond cleavage allowed for each
%             fragmentation mode.
%         nmFrag: max. number of non-glycan PTM detachment allowed for each
%             fragmentation mode.
%         maxLag: max. lag for cross correlation calculation.
%         CutOffMed: denoising power, 1 of 2. See SPECTRADENOISING for
%             details.
%         FracMax: denoising power, 2 of 2. See SPECTRADENOISING for
%             details.
%         selectpeak: specific m/z values to search in each MS2 scan.
% msdatatype: result file type, mzXML or DTA or .mat.
% datatyp: the data type of each field in scoredata, "numeric" (data is a
%     numeric array), "cell-numeric" (data is cell array of numebrs), "string" (data
%     is cell array of strings).
% scoredata: structure, scores in the input file. The field names are
%     direct copy of column headers. If the header is not in the list below
%     (case insensitive), the data will not appear in this output.
%     Acceptable header names are:
%         (numeric): "SCAN", "EXPT", "MONO", "THEO", "MDIFF", "CHARGE",
%             "PEAKLAG", "HTCENTER", "HTAVG", "PERCENTIONMATCH", "PVALUE",
%             "DECOYRATIO", "TOP10", "NMFRAG", "NGFRAG", "NPFRAG", "ENSCORE",
%             "DECOYENSCORE", "QUANT", "RETIME", "GLYCOV", "FRACPEPFRAGMATCH",
%             "FRACGLYFRAGMATCH"
%         (cell-numeric): "SELECTPEAK", "PROTEINID", "Y0Y1Y2"
%         (string): "PROTEIN", "SGP", "FRAGMODE"
% 
% NOTE:
% N/A
% 
% Example:
% [headerpar,msdatatype,datatyp,scoredata]=scoreCSVread('test_score.csv');
% 
% Children function: 
% READTABULARDATA  READPARFROMHEADER  
% See Also: 
% SCORECSVWRITE

% Author: Gang Liu
% Date Lastly Updated: 11/10/20 by Kai Cheng

[scoredata,rawdatabyline,linenum] = readtabulardata(fname);

% The section where the header is stored. Reason for minus one is because
%     "linenum" includes the column names of the score table.
headerpart = rawdatabyline(1:linenum - 1);  
[headerpar,msdatatype] = readparfromheader(headerpart);

% Know the data type of each column of the score table
varisnumeric = ResultItem.itemisnumeric;  
varisnumcell = ResultItem.itemisnumcell;
varisstring = ResultItem.itemisstring;

% Retrieve the column names, so different data types can be treated
%     accordingly
fldnms = fieldnames(scoredata);
datatyp = cell(1,length(fldnms));
keepind = ~cellfun(@isempty,scoredata.(fldnms{1}));
for i = 1:length(fldnms)
    scoredata.(fldnms{i}) = scoredata.(fldnms{i})(keepind);
end
for i = 1:length(fldnms)
    
    % Numeric types are saved directly
    if ismember(upper(fldnms{i}),varisnumeric)
        scoredata.(fldnms{i}) = cellfun(@str2double,scoredata.(fldnms{i}));
        datatyp{i} = 'numeric';
        
    % Cell-numeric types are texts originally, convert them to numbers then
    %     saved in cell array
    elseif ismember(upper(fldnms{i}),varisnumcell)
        scoredata.(fldnms{i}) = cellfun(@str2num,scoredata.(fldnms{i}),'UniformOutput',false);
        datatyp{i} = 'cell-numeric';
        
    % Cell array of strings are saved directly
    elseif ismember(upper(fldnms{i}),varisstring)
        tempdata = scoredata.(fldnms{i});
        for j = 1:length(tempdata)
            if strcmp(tempdata{j}(1),'"') && strcmp(tempdata{j}(end),'"')
                tempdata{j} = tempdata{j}(2:end-1);
            end
        end
        scoredata.(fldnms{i}) = tempdata;
        datatyp{i} = 'string';
    end
end
end

function [headerpar,msdatatype] = readparfromheader(headercell)
% READPARFROMHEADER: read header information and return a structure,
%     field names are name of information, values are information
%     themselves
% Input:
% headercell: n x 1 cell array of strings. The contents of the .csv file
%     that contains header information
% 
% Output:
% headerpar: structure, header information names and information themselves
% msdatatype: string, the type of experiment data.

if ~isempty(headercell)
    for ii = 1:1:length(headercell)
        
        % Due to the fact that .csv files may be modified by text
        %     processing softwares (e.g. Excel), which may add extra commas, 
        %     these extra commas must be recognized and removed before
        %     their contents can be read correctly
        if strcmpi(headercell{ii,1},',')
            headercell{ii,1} = headercell{ii,1}(2:end);
        end
        CommaPos=regexp(headercell{ii},',');
        if(length(CommaPos)>=2)
            itemname = headercell{ii}(1:CommaPos(1) - 1);
        elseif isempty(CommaPos)
            itemname = headercell{ii};
        else
            itemname = headercell{ii}(1:CommaPos(1) - 1);
        end
        CommaPos=regexp(headercell{ii},',');
        if(length(CommaPos)>=2)
            itemcontent = headercell{ii}(1:CommaPos(1) - 1);
        elseif isempty(CommaPos)
            itemcontent = headercell{ii};
        else
            itemcontent = headercell{ii}(CommaPos(1) + 1:end);
        end
        itemname = strrep(itemname,' ','_');
        itemname = erase(itemname,{':','.','=','(',')'});
        headerpar.(itemname) = itemcontent;
    end
    
    % Specific field names that are important 
    if isfield(headerpar,'ExperimentData')
        
        % Full path and name of experiment data file
        exptdatafilename = headerpar.ExperimentData;
    elseif isfield(headerpar,'mat_File_Name')
        exptdatafilename = headerpar.mat_File_Name;
    end
   [~,~,ext] = fileparts(exptdatafilename);
   
   % Experiment data file type
   switch upper(ext)
       case '.MAT'
           msdatatype = 'mat';
       case '.MZXML'
           msdatatype = 'mzml';
       case '.MZML'
           msdatatype = 'mzxml';
   end
else
    % Sometimes user removes the header, in which case empty values will be
    %     returned
    headerpar = [];
    msdatatype = '';
end
end

function [csvstruct,databyline,line] = readtabulardata(fname)
% READTABULARDATA: read the score part of the file, return structure
%     containing score table. Structure field names are direct copy of 
%     header, data are stored in cell array of strings. 
% 
% Output:
% csvstruct: structure with score data.
% databyline: cell array of strings, the whole content of the file, line by
%     line.
% line: the row number of the score table header, therefore file header is
%     stored in row 1 through row (line - 1)
% 
% Note: 
% This function decides where the scores start by finding keywords
%     "Enscore", "AUC", "Scan" or "Charge", if any one of the above
%     appears, the row where it was found will be treated as the start of
%     the score table

% Load the file, show error information when file is not available
fid = fopen(fname,'r');
if fid == -1
    error('Cannot load file');
end

% Read the entire file row by row
databyline = textscan(fid,'%s','delimiter','\n');
fclose(fid);

% Using the argument above, result from TEXTSCAN must be unpacked first
databyline = databyline{1};  
line = 1;

% Search row by row, looking for key words indicating score table column
% names. It stops at that row
while ~(~isempty(strfind(databyline{line},'Enscore')) || ...
        ~isempty(strfind(databyline{line},'AUC')) || ...
        ~isempty(strfind(databyline{line},'Scan')) || ...
        ~isempty(strfind(databyline{line},'Charge')))
    line = line + 1;
end

% In .csv files, data are separated by commas
columnames = strsplit(databyline{line},',');
columnames = strrep(columnames,'"','');
columnames = columnames(~cellfun(@isempty,columnames));

% The data in the columns mentioned below are in text format. Others are
%     treated as numbers.
specialfields = {'Protein','%s';...
    'SGP','%s';...
    'Fragmode','%s';...
    'SelectPeak','%s'};
fieldformat = repmat({'%s'},1,length(columnames));
for i = 1:length(columnames)
    isspecialfldnm = ismember(upper(specialfields(:,1)),upper(columnames{i}));
    if any(isspecialfldnm)
        fieldformat{i} = specialfields{isspecialfldnm};
    end
end

% Load the score table
datatable = cell(size(databyline,1) - line,length(columnames));
for i = 1:size(databyline,1)-line
    thislinetext = databyline{line + i};
    [optcontent,optcontentstart,optcontentend] = regexp(thislinetext,'["''“].*?["''”]','match','start','end');
    tempenclosedseq = thislinetext;
    if ~isempty(optcontent)  % block customization info ("noise")
        for j = 1:length(optcontentstart)
            tempenclosedseq(optcontentstart(j):optcontentend(j)) = blanks(optcontentend(j) - optcontentstart(j) + 1);
        end
    end
    if length(strfind(tempenclosedseq,',')) < length(columnames)
        tempenclosedseq = [tempenclosedseq,','];
    end
    [thisline,commapos] = regexp(tempenclosedseq,'.*?,','match','end');
    commapos = commapos - 1;
    for j = 1:length(thisline)
        if any(commapos(j) == optcontentend)
            thisline{j} = optcontent{optcontentend == commapos(j)}(2:end-1);
        else
            thisline{j} = thisline{j}(1:end-1);
        end
    end    
    for j = 1:length(columnames)
        temp = char(thisline(j));
        if strcmp(fieldformat{j},'%s')
            if length(temp) > 2 && strcmp(temp(1),'"') && strcmp(temp(end),'"')
                temp = temp(2:end-1);
            end
        end
        datatable{i,j} = temp;
    end
end
for i = 1:size(datatable,2)
    if strcmp(fieldformat{i},'%f')
        datatable(:,i) = cellfun(@str2double,datatable(:,i),'UniformOutput',false);
    end
end
for i = 1 : length(columnames)
    csvstruct.(columnames{i}) = datatable(:,i);
end
end