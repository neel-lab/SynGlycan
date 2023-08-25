function fixedptmopt = fixedptmread(fptmfullfilename)
%fixedptmread: Load the option for fixed peptide post-translational
%      modification.
%
% Syntax:
% fixedptmopt = fixedptmread(fptmfullfilename)
%
% Input:
% fptmfullfilename: the full file name including the path.
%
% Output:
% fixedptmopt: a 1 x n structure with the following fields:
%     aaresidue: string. The amino acid to be modified.
%     new: string. The product of the modification.
%     protpos: 1 x n numerical array. The position of the amino acid to be
%         modified.
%     mod: string. The modification itself (the content inside "< >")
%
% Note:
% For a detailed description of fixed PTM file format, please read the
%     example file in "Test Data Folder'\Protein Digestion\FPTM.txt"
%
% Examples:
% fixedptmopt = fixedptmread('"Test Data Folder"\Protein Digestion\FPTM.txt')
%
% fixedptmopt =
%
%   struct with fields:
%
%     aaresidue: 'C'
%           new: 'C<i>'
%       protpos: 0
%           mod: 'i'
%
% See also:
% FASTAREAD  VARPTMREAD  CLEAVEPROT  DIGESTSGP

% Author: Gang Liu
% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
%(c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 03/01/19

fid = fopen(fptmfullfilename);  % Open fixed PTM file
if(fid == -1)  % File not available
    fclose(fid);
    error('MATLAB:GLYCOPAT:FILENOTFOUND','file is not found in the search path');
end
alltext = textscan(fid,'%s','Delimiter',{'\n';'\r';'\r\n';'\n\r'});  % Read file content
fclose(fid);  % Release the file after reading is complete
alltext = alltext{1};  % Function "textscan" wraps the result in a cell, now unfold it
alltext = alltext(~cellfun(@isempty,alltext));  % Remove empty lines
alltext = cellfun(@strtrim,alltext,'Uniformoutput',false);  % Remove leading and trailing whitespaces of each line

% Remove comment lines: comment lines start with a percentage mark "%"
iscomment = zeros(size(alltext));  
for Ii = 1:length(alltext)
    if strcmp(alltext{Ii}(1),'%')
        iscomment(Ii) = 1;
    end
end
alltext_nocomment = alltext(~logical(iscomment));
alltext_nocomment = alltext_nocomment(~cellfun(@isempty,alltext_nocomment));

% Standard fixed PTM format: [Amino acid]<Modification>, e.g. "C<i>"
pat = '(?<aa>[a-z_A-Z])<(?<mod>[a-z_A-Z0-9]+)>';  

% Read each line. Starting from last one can avoid format errors.
for Ii = length(alltext_nocomment):-1:1
    splittedline = strsplit(alltext_nocomment{Ii},' ');  % The 3 parts describing each PTM is separated by whitespaces
    tempprotpos = str2num(splittedline{2});  % Position of amino acid to be modified
    fixoptnames = regexp(splittedline{3},pat,'names');  % Check the 3rd part
    % The product of fixed PTM may contain only 1 amino acid, e.g. "C<i>" is an acceptable format,
    %     but "MC<i>" is not allowed. These amino acids must be in the
    %     predefined list as described in file "Aminoacid.m"
    
    % Error detection - more than 1 AA or AA is not listed in AMINOACID
    if (length(fixoptnames.aa) ~= 1) || ~Aminoacid.isaastring(fixoptnames.aa)
        error('MATLAB:GLYCOPAT:NOTVALIDAACHAR','NOT VALID AMINO ACID CHARACTER');
    end
    
    % Check position of amino acid to be modified
    if length(tempprotpos) > 1  % if number of position > 1, none shall be "0"
        if any(tempprotpos == 0)
            error(['Fixed modification position information error in: "' alltext_nocomment{Ii} '"']);
        end
    elseif length(tempprotpos) == 1 && tempprotpos == 0  % if number is "0", the amino acid letter
        % shall not be "X"
        if strcmpi(splittedline{1},'X') && ~ismember(fixoptnames.mod,Modification.isotagname)
            error('Cannot set "X" at position "0" in fixed modification')
        end
    end
    fixedptmopt(Ii).aaresidue = splittedline{1};
    fixedptmopt(Ii).new = splittedline{3};
    fixedptmopt(Ii).protpos = tempprotpos;
    fixedptmopt(Ii).mod = fixoptnames.mod;
end
if isempty(alltext_nocomment)
    fixedptmopt.aaresidue = 'X';
    fixedptmopt.new = 'X';
    fixedptmopt.protpos = 0;
    fixedptmopt.mod = ' ';
end
end