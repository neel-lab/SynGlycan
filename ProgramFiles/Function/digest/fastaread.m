function protdata = fastaread(protfullfilename)
% FASTAREAD: Read protein sequence from a sequence file in FASTA format
%
% Syntax:
% protdata = fastaread(protfullfilename)
%
% Input:
% pepfilefullname: the full path for a protein sequence file.
%
% Output:
% protdata: 1 x n structure with fields: sequence, header.
%     sequence: string, amino acid sequence of protein;
%     header: string, FASTA header of protein.
%
% Examples:
% protdata = fastaread('"Test Data Folder"\Protein Digestion\10000PROT.fasta')
%
% protdata =
%
%   10000×1 struct array with fields:
%
%     header
%     sequence
%
% See also:
% FIXEDPTMREAD VARPTMREAD CLEAVEPROT DIGESTSGP

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
%(c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 03/01/19

% Try to open the file, assign a handle for content access.
fid = fopen(protfullfilename);
if(fid==-1) % File not found
    error('MATLAB:GlycoPAT:FILENOTFOUND','File is not found in specified path.');
end

% Decide file type through file extension.
% If sequence file type is supported, start reading the file.
[~,~,fileformat] = fileparts(protfullfilename);
if(isequal(fileformat,'.txt') || isequal(fileformat,'.fasta') ...
        || isequal(fileformat,'.fas') || isequal(fileformat,'.fa') ...
        || isequal(fileformat,'.seq') || isequal(fileformat,'.fsa') ...
        || isequal(fileformat,'.faa'))
    
    % Read sequence from file.
    try
        ftext = textscan(fid,'%s','delimiter','\n');  % Read all text in the file.
        fclose(fid);  % Close file after data has been read.
        ftext = ftext{:};  % Rearrange data for easier processing.
        ftext = ftext(~cellfun(@isempty,ftext));  % Remove empty rows in file.
        ftext = cellfun(@strtrim,ftext,'Uniformoutput',false);  % Remove leading and trailing whitespaces.
        iscommentline = false(size(ftext));
        % Assign a "type" to each line, for distinguishing comments
        %     from real data. Lines started with "%" are comments, they will not be included as
        %     part of the analysis.
        for i = 1:length(ftext)
            if strcmpi(ftext{i}(1),'%')
                iscommentline(i) = true;
            end
        end
        ftext = ftext(~iscommentline);  % Remove comment lines.
    catch readingErr  % If reading goes wrong ("readingErr" is a built-in error type of MATLAB).
        if strcmpi(readingErr.identifier,'MATLAB:nomem')
            error('MATLAB:GlycoPAT:NOMEM','FILE IS TOO BIG');
        else  % Other than memory overflow, display the whole error message to user.
            rethrow(readingErr);
        end
    end
    
    % In FASTA format protein names are marked by ">" at the beginning of
    %     the line. This is the sole indicator of protein names
    %     recognizable by GlycoPAT.
    headerLines = strncmp(ftext,'>',1);    
    if ~any(headerLines)  % No protein name found.
        error('MATLAB:GlycoPAT:FASTAREADERROR','NO HEADER LINES IN FASTA FILE');
    end
    
    numSeqs = sum(headerLines);  % Number of ">" equals to the number of headers, 
    %    which equals to the number of proteins.
    seqStarts = [find(headerLines);size(ftext,1) + 1];
    % The latter part ("size(ftext,1) + 1") is a "pseudo header". Anything between the two
    %     headers is the protein amino acid sequence. This "pseudo header" allows
    %     the last protein be read.
    
    % Note: If file coding is Unicode, at the beginning of the file there
    %     will be a byte order mark (BOM). BOM will scramble the protein
    %     name recognization. In that case, the code below will take
    %     whatever in the first line as the protein name.
    if seqStarts(1) ~= 1
        seqStarts = [1;seqStarts];
        numSeqs = numSeqs + 1;
    end
    protdata(numSeqs,1).header = '';
    
    unkprotind = 1;
    try
        for theSeq = 1:numSeqs
            thesymbol = strfind(ftext{seqStarts(theSeq)},'>');  % Check for > symbol
            protdata(theSeq).header = ftext{seqStarts(theSeq)}(thesymbol(1) + 1:end);
            % convert 1x0 empty char array to '';
            if isempty(protdata(theSeq).header)  % Protein name is empty but does have ">"
                protdata(theSeq).header = ['Incognito Protein ',num2str(unkprotind)];
                unkprotind = unkprotind + 1;
            end
            firstRow = seqStarts(theSeq) + 1;
            lastRow = seqStarts(theSeq + 1) - 1;  % Any text between two protein headers
            numChars = cellfun('length',ftext(firstRow:lastRow));  % In case sequence consists of multiple lines of text
            sequence = blanks(sum(numChars));  % Apply for space in memory first, then fill with text.
            pos = 1;
            for i=firstRow:lastRow
                str = ftext{i};
                len =  length(str);
                if len == 0
                    break
                end
                sequence(pos:pos+len-1) = str;
                pos = pos+len;
            end
            sequence = upper(strtrim(sequence));  % Amino acids letters are uppercase
            % PEPTIDE SEQUENCE CHARACTER CHECK
            aaexpr = ['[^' Aminoacid.getaachar ']'];
            anynonletterchar = regexp(sequence,aaexpr,'once');
            if(~isempty(anynonletterchar))
                error('MATLAB:GlycoPAT:NONPEPSEQ',['INCORRECT PEPTIDE SEQUENCE CHARACTER: "' ...
                    sequence(anynonletterchar) '" IN PROTEIN: "' protdata(theSeq).header '"']);
            end
            protdata(theSeq).sequence = sequence;
        end
    catch allExceptions
        error(message('MATLAB:GlycoPAT:INCORRECTFASTAFORMAT'))
    end
else
    error('MATLAB:GlycoPAT:NONSUPPORTEDFILETYPE','FILE TYPE IS NOT SUPPORTED');
end
end