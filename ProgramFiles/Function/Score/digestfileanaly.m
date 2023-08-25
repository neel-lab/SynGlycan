function allpepdata = digestfileanaly(Pepfile,pepfileversion)
% DIGESTFILEANALY: Analysis digested glycopeptide file to present its
% content
%
% Syntax:
% allpepdata = digestfileanaly(Pepfile,pepfileversion)
%
% Input:
% Pepfile: string, location of digested peptide file.
% pepfileversion: number, currently "1" is supported by GlycoPAT2
%
% Output:
% allpepdata: structure, name of the fields is associated with contents in
%     the file. Available names and their corresponding in file identifiers are:
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%     Field name                        Format                                 Description
%     FASTAhead                       n x 1 cell array of              The header part of each protein in the FASTA file.
%                                                      strings
%     glypep                                n x 1 cell array of              Digested glycopeptides, each n is a digested protein, each m is a
%                                                      m x 1 cell array of             peptide backbone ofglycopeptide. See FINDMS1MATCH for
%                                                      strings                                  detail.
%                                                      
%     ptminfo                              Structure                            Contains 2 fields: "mod" and "uniind". "mod" is a n x 1 cell array
%                                                                                                    of strings, containing all PTMs that may appear on glycopeptides.
%                                                                                                    "uniind" is n x 1 double, the locator of each PTM. See
%                                                                                                    ASSEMBLEGP for the explaination of locator.
%     (Below are not essential to scoring, they are here for record keeping)
%     summary                           n x 1 cell array of              The summary of protein digestion operation.
%                                                      strings
%     enzyme                              n x 1 cell array of              The enzyme(s) used in protein digestion.
%                                                      strings
%     fixedptm                            Structure                             Fixed PTM information.
%         ~.fixptmuniind             n x 1 double                        Locator of fixed PTM.
%         ~.aaresidue                   n x 1 cell array of               Modified amino acid residue.
%                                                      strings
%         ~.mod                             n x 1 cell array of               Sequence of fixed PTM.
%                                                      strings
%         ~.peppos                        n x 1 double                        Position of fixed PTM on peptide.
%     varptm                               Structure                             Variable PTM information.
%         ~.aaresidue                   n x 1 cell array of               Modified amino acid residue.
%                                                      strings
%         ~.mod                             n x 1 cell array of               Sequence of variable PTM.
%                                                      strings
%         ~.pos                               n x 1 cell array of               Position of variable PTM on peptide.
%                                                       double
%         ~.varptmuniind            n x 1 double                        Locator of variable PTM.
%     minptm                              Double                                 Minimum number of PTMs allowed on each digested glycopeptide.
%     maxptm                             Double                                 Maximum number of PTMs allowed on each digested glycopeptide.
%     minpeplen                         Double                                 Minimum length of digested glycopeptide.
%     maxpeplen                        Double                                 Maximum length of digested glycopeptide.
%     missedmax                        Double                                 Maximum number of missed cleavage allowed.
%     nocleave                            n x 1 cell array of               The FASTA header of proteins that cannot be digested under
%                                                      strings                                  present settings.
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%
% Note:
% A typical digested peptide file have multiple sections. Each section has
%     a header marked by "#" at the beginning. Currently available sections
%     are:
%         #SUMMARY: number of proteins digested, number of peptides
%             generated, time consumed.
%         #VERSION: version of this file's formatting.
%         #PROTSEQPATH: the full file name and path of the .fasta file
%             used.
%         #ENZYME: proteases used.
%         #PTM_FIXED: fixed PTMs used.
%         #PTM_VARIABLE: variable PTMs used.
%         #PTM_MIN: minimum number of PTMs allowed on peptide.
%         #PTM_MAX: maximum number of PTMs allowed on peptide.
%         #PEPLEN_MIN: minimum length of peptide allowed.
%         #PEPLEN_MAX: maximum length of peptide allowed.
%         #MISSEDMAX: maximum number of missed cleavages allowed.
%         #PROTEIN_NOT_DIGESTED: name of the proteins that generated no
%             peptide.
%         #SEQUENCES: sequence of digested peptides. Peptides are grouped
%             by proteins, each protein starts with its FASTA header
%             starting with ">".
%         #END: end of the file.
%
% Example:
% N/A
%
% Children function:
% N/A

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 04/01/2020

fileID = fopen(Pepfile,'r');  % Open the file with permission "Read Only"
if fileID == -1
    error('DIGEST FILE CANNOT BE READ');
end
texts = textscan(fileID,'%s','Delimiter',{'\n';'\r';'\r\n';'\n\r'});  % Read the file row by row
fclose(fileID);  % Release the file when all data has been read
texts = texts{1};  % TEXTSCAN packs file content in a 1 x 1 cell, here it's unpacked
texts = texts(~cellfun(@isempty,texts));  % Remove empty lines
texts = cellfun(@strtrim,texts,'UniformOutput',false);  % Remove leading and trailing whitespaces
iscomment = zeros(size(texts));  % Lines that are comments start with "%" 
for i = 1:length(texts)
    if strcmp(texts{i}(1),'%')  % Line starts with "%"
        iscomment(i) = 1;
    end
end
alltext_nocomment = texts(~logical(iscomment));  % Comment free text
texts = alltext_nocomment(~cellfun(@isempty,alltext_nocomment));  % Remove empty lines - doublecheck
isheader = zeros(size(texts));  % Digestion file is divided into several sections. Each section has a header, started by "#"
for i = 1:length(texts)  % Find the lines that are headers
    if strcmp(texts{i}(1),'#')
        isheader(i) = 1;
    end
end
sectionstartend = zeros(sum(isheader)-1,2);  % Section starting and ending row number
isheader = find(isheader);
sectionstartend(:,1) = isheader(1:length(isheader)-1) + 1;
sectionstartend(:,2) = isheader(2:length(isheader)) - 1;  % [Section_1_start, Section_1_end;Section_2_start, Section_2_end;...]
switch pepfileversion
    case 1  % GlycoPAT 2.0 format. 
        % All modification sites on glycopeptides that are marked will carry a PTM. 
        % For glycopeptide carrying variable PTMs, it is possible that
        %     mod. site be vacant, in which case 2 glycopeptides are needed:
        %     one for vacant, the other for occupied
        for i = 1:length(isheader) - 1  % Minus 1 because the final section has another "header" called "#END"
            data = texts(sectionstartend(i,1):sectionstartend(i,2));
            switch texts{isheader(i)}
                case '#SUMMARY'
                    allpepdata.summary = data;  % Digestion parameters: min. max. peptide length,...
                case '#VERSION'
                    allpepdata.version = data;  % Digestion parameters: min. max. peptide length,...
                case '#PROTSEQPATH'
                    allpepdata.protseqpath = data;  % Digestion parameters: min. max. peptide length,...
                case '#ENZYME'
                    allpepdata.enzyme = data;  % Protease used and regular expression for cleavage behavior
                case '#PTM_FIXED'
                    fixedptmdat = cell(size(data,1),4);  % List of fixed PTM
                    for j = 1:size(data,1)
                        fixedptmdat(j,:) = strsplit(data{j},' ');  % Format: {Marker} {Original AA} {Modified AA} {Modification sites}
                        % For example: {3} {C} {C<i>} {0} means Cystine
                        %     carboxymethylation will appear as C{3} in
                        %     digestion library. Its final form will be C<i>
                    end
                    allpepdata.fixedptm.fixptmuniind =  cellfun(@str2double,fixedptmdat(:,1));  % Originally "Marker" is in text form, is converted to number here
                    allpepdata.fixedptm.aaresidue = fixedptmdat(:,2);
                    allpepdata.fixedptm.mod = cellfun(@setdiff,fixedptmdat(:,3),fixedptmdat(:,2),...
                        repmat({'stable'},size(fixedptmdat,1),1),'uniformoutput',false);  % Get the modification itself, e.g. <i> from C<i>
                    allpepdata.fixedptm.peppos =  cellfun(@str2double,fixedptmdat(:,4));  % Fixed PTM modification position
                    allpepdata.fixedptm.isglycan =  ~cellfun(@isempty,strfind(allpepdata.fixedptm.mod,'{'));
                case '#PTM_VARIABLE'
                    varptmdat = cell(size(data,1),4);  % Format is identical to fixed PTM
                    for j = 1:size(data,1)
                        varptmdat(j,:) = strsplit(data{j},' ');
                    end
                    allpepdata.varptm.aaresidue = varptmdat(:,2);
                    allpepdata.varptm.mod = varptmdat(:,3);
                    allpepdata.varptm.pos = varptmdat(:,4);
                    allpepdata.varptm.varptmuniind = cellfun(@str2double,varptmdat(:,1));
                    allpepdata.varptm.isglycan =  ~cellfun(@isempty,strfind(allpepdata.varptm.mod,'{'));
                % See documentation for the meaning of these sections below
                case '#PTM_MIN'
                    allpepdata.minptm = str2double(data{1});
                case '#PTM_MAX'
                    allpepdata.maxptm = str2double(data{1});
                case '#PEPLEN_MIN'
                    allpepdata.minpeplen = str2double(data{1});
                case '#PEPLEN_MAX'
                    allpepdata.maxpeplen = str2double(data{1});
                case '#MISSEDMAX'
                    allpepdata.missedmax = str2double(data{1});
                case '#PROTEIN_NOT_DIGESTED'
                    allpepdata.nocleave = data;
                    
                % Reading glycopeptides
                case '#SEQUENCES'
                    isheaderline = zeros(size(data));  % here "header" is the protein name
                    for j = 1:length(data)
                        if strcmp(data{j}(1),'>')  % FASTA format requires protein names start with ">"
                            isheaderline(j) = 1;
                        end
                    end
                    headerlinepos = find(isheaderline);
                    FASTAhead = cell(size(headerlinepos));  % FASTAhead stores protein names in FASTA file
                    for j = 1:length(headerlinepos)
                        FASTAhead{j} = data{headerlinepos(j)}(2:end);
                    end
                    allpepdata.FASTAhead = FASTAhead;
                    sectionind = [headerlinepos(1:end) + 1,[headerlinepos(2:end)-1;length(data)]];  % Starting and finishing row numbers for glycopeptides
                    glypep = cell(size(sectionind,1),1);
                    for j = 1:size(sectionind,1)
                        glypep{j} = data(sectionind(j,1):sectionind(j,2));
                    end
                    allpepdata.glypep = glypep;
            end
        end
        
        % Store variable and fixed PTM into the same variable.
        % During digestion, these two types of PTM were given unique codes,
        %     their modification sites were already embedded in
        %     glycopeptide sequences, which satifies rules of modification.
        %     It is OK to put them into one variable.
        ptminfo.mod = {};
        ptminfo.uniind = [];
        ptminfo.isglycan = [];
        ptminfo.isvarmod = [];
        if isfield(allpepdata,'varptm')
            ptminfo.mod = allpepdata.varptm.mod;
            ptminfo.uniind = allpepdata.varptm.varptmuniind;
            ptminfo.isglycan = allpepdata.varptm.isglycan;
            ptminfo.isvarmod = true(size(ptminfo.mod));
        end
        if isfield(allpepdata,'fixedptm')
            ptminfo.mod = [ptminfo.mod;allpepdata.fixedptm.mod];
            ptminfo.uniind = [ptminfo.uniind;allpepdata.fixedptm.fixptmuniind];
            ptminfo.isglycan = [ptminfo.isglycan;allpepdata.fixedptm.isglycan];
            ptminfo.isvarmod = [ptminfo.isvarmod;false(size(allpepdata.fixedptm.mod))];
        end
        allpepdata.ptminfo = ptminfo;
    case 2
        % Proposed new format: Glycopeptide will have all potential
        %     modification sites marked. If PTM is variable, it may not be
        %     occupied. This reduces file size. 
        %     NOT AVAILABLE IN CURRENT VERSION
end
end