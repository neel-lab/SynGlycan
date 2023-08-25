function varptmopt = varptmread(vptmfullfilename)
%VARPTMREAD: Load the option for variable peptide post-translational
%    modification.
%
% Syntax:
% vptmopt = varptmread(varptmfilefullname)
%
% Input:
% vptmfullfilename: the full file name including the path.
%
% Output:
% varptmopt: a 1 x n structure with the following fields:
%     mod: string. The modification itself.
%     aaresidue: string. The amino acid to be modified.
%     protpos: number. The position of the amino acid to be modified.
%     isglymod: logical. Whether the modification is a glycan.
%
% Note:
% For a detailed description of fixed PTM file format, please read the
%     example file in "Test Data Folder'\Protein Digestion\VPTM_N79+NG1.txt"
%
% Examples:
% varptmopt = varptmread('"Test Data Folder"\Protein Digestion\VPTM_N79+NG1.txt')
%
% varptmopt =
%
%   1×80 struct array with fields:
%
%     mod
%     aaresidue
%     protpos
%     isglymod
%
% See also:
% FASTAREAD FIXEDPTMREAD CLEAVEPROT DIGESTSGP MODIFICATION

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 04/01/2020

fid = fopen(vptmfullfilename);  % Open the file for reading
if fid == -1  % if file not found
    fclose(fid);
    error('MATLAB:GLYCOPAT:FILENOTFOUND','file is not found in the search path');
end

alltext = textscan(fid,'%s','Delimiter',{'\n';'\r';'\r\n';'\n\r'});  % Read the file
fclose(fid);  % Release the file after reading is complete
alltext = alltext{1};  % Function "textscan" wraps the result in a cell, now extracting it
alltext = alltext(~cellfun(@isempty,alltext));  % Remove empty lines
alltext = cellfun(@strtrim,alltext,'UniformOutput',false);  % Remove leading and trailing whitespaces
iscomment = zeros(size(alltext));  % Remove comment lines: comment lines start with a percentage mark "%"
for i = 1:length(alltext)
    if strcmp(alltext{i}(1),'%')
        iscomment(i) = 1;
    end
end
alltext_nocomment = alltext(~logical(iscomment));
modif = cell(size(alltext_nocomment));  % prepare containers for information
aaresidue = modif;
pos = modif;
for i = 1:length(alltext_nocomment)
    splittedline = strsplit(alltext_nocomment{i},' ');  % The 3 parts describing each PTM is separated by whitespaces
    modif{i} = splittedline{1};  % Modification itself
    aaresidue{i} = splittedline{2};  % Amino acid to be modified
    pos{i} = str2num(splittedline{3});  % Position of amino acid to be modified
end
varptmopt = struct('mod',[],'aaresidue',[],'protpos',[],'isglymod',[]);

% Non-glycan PTM regular expression
pat = '<[a-z_A-Z]>';
pat_name = '<(?<mod>[a-z_A-Z])>';

nummod = length(aaresidue);
mod1let = Modification.mod1let;  % Retrieve single letter modification codes
mod1let = [mod1let{:}];
nonmod1let = ['[^',mod1let,']'];
for i = nummod:-1:1
    modifstring = modif{i};
    modstringfound = regexp(modifstring,pat,'once');  % Variable is not empty if
    % it is non-glycan PTM, represented by "< >"
    isglycanmod = isempty(modstringfound);
    if ~isglycanmod  % Non-glycan PTM
        varoptmod = regexp(modifstring,pat_name,'names');
        % Check modification letter. Only pre-defined single letter code
        %     can be used
        if length(varoptmod.mod) == 1  % Check whether it is single letter
            if(~isempty(regexpi(varoptmod.mod,nonmod1let,'once')))  % Check whether it is pre-defined
                error('MATLAB:GlycoPAT:INVALIDMODCHAR','INVALID MODIFICATION LETTER');
            end
        else  % Multiple letter as code - ERROR
            error('MATLAB:GlycoPAT:INVALIDMODCHAR','INVALID MODIFICATION LETTER');
        end
        varptmopt(i).mod = ['<' varoptmod.mod '>'];
        varptmopt(i).isglymod = false;
    else
        varptmopt(i).mod = modifstring;
        varptmopt(i).isglymod = true;
    end
    if strcmpi(aaresidue{i},'X') && pos{i,1} == 0
        error('Invalid amino acid and position combination: " X 0"')
    end
    varptmopt(i).aaresidue = aaresidue{i};
    varptmopt(i).protpos = pos{i,1};  % might be multiple numbers so have to use str2num
end

for i = 1:nummod
    % Check amino acid letter
    if length(varptmopt(i).aaresidue) < 1 || ~isempty(regexpi(varptmopt(i).aaresidue,...
            '[^ARNDCQEGHILKMFPSTWYVBZX,]'))
        error('MATLAB:GlycoPAT:NOTVALIDAACHAR','NOT VALID AMINO ACID CHARACTER');
    end
    % Check position
    if ~isnumeric(varptmopt(i).protpos) || isempty(varptmopt(i).protpos) || sum(isnan(varptmopt(i).protpos)) >= 1
        error('MATLAB:GlycoPAT:NOTVALIDPOS','NOT VALID MODIFICATION POSITION');
    end
end
end