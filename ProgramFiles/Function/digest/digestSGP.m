function [finalallprod,outputstring,reportstring]=digestSGP(protseq,protseqheader,enzname,digestopt)
% DIGESTSGP: Create a list of enzyme digested fragments given the protein
%     sequence, list of glycan and non-glycan modifications, and other
%     digestion parameters.
%
% Syntax:
% [finalallprod,outputstring,reportstring]=digestSGP(protseq,protseqheader,enzname,digestopt)
%
% Input:
% protseq: n x 1 cell array of strings. FASTA sequence of proteins to be
%     digested
% protseqheader: n x 1 cell array of strings. Name of proteins
% enzname: 1 x n cell array of strings. Name of the protease(s)
% digestopt: structure. Digestion options. Fields are:
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%     Field name                        Format                                 Description
%     missedmax                        Double                                 Maximum number of missed cleavages allowed.
%     minpeplen                        Double                                  Minimum length of digested peptide.
%     maxpeplen                        Double                                 Maximum length of digested peptide.
%     minptm                              Double                                 Minimum number of variable PTMs allowed on peptide.
%     maxptm                             Double                                 Maximum number of variable PTMs allowed on peptide.
%     fixedptm                            1 x n structure                   List of variable PTM. See FIXEDPTMREAD for detail.
%     varptm                               1 x n structure                   List of variable PTM. See VARPTMREAD for detail.
%     dispprog                            Logical                                  Whether or not display a wait bar to show progress.
%     protfilepath                      String                                    Full file path of the .fasta file containing protein sequences.
%     doparacomp                    Double                                  Number of CPU cores to use.
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%
% Output:
% finalallprod: structure. Digestion results. See children function
%     RESULT2TEXT for details.
%     
% outputstring: string. Formatted digestion result ready to be written to
%     text file. See children function RESULT2TEXT for detail.
% reportstring: n x 1 cell array of strings. The summary of program's
%     result. The contents are: number of  proteins digested; number of
%     generated peptides; which proteins created no peptides after
%     digestion and time consumed.
%
% Note:
% Maximum number of modifications allowed is 65429.
%
% Examples:
% N/A. Set breakpoints in DIGESTGUI, load the test files, then run.
%
% Children Function:
% CLEAVEPROTEINS ADDVARTOPEP FORMATVARIABLE
%     APPLYVARMOD PERMUTATION RESULT2TEXT
%
% See also:
% DIGESTGUI  FASTAREAD  FIXEDPTMREAD  VARPTMREAD  CLEAVEPROT

% Author: Sriram Neelamegham, Gang Liu
% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 04/01/2020

tic
finalallprod=[];
dispprogress = digestopt.dispprog;
%% HANDLING INPUT ARGUMENT
% Check whether options for digestion are acceptable
if ( ((~isempty(digestopt.missedmax)) && (~isnumeric(digestopt.missedmax)))) || ...
        ((~isempty(digestopt.minpeplen)) && (~isnumeric(digestopt.minpeplen))) || ...
        ((~isempty(digestopt.maxpeplen)) &&  (~isnumeric(digestopt.maxpeplen))) || ...
        ((~isempty(digestopt.minptm)) && (~isnumeric(digestopt.minptm)))  || ....
        ((~isempty(digestopt.maxptm)) && (~isnumeric(digestopt.maxptm)))
    error('MATLAB:GlycoPAT:INCORRECTTYPE','Wrong input type for the propertyname');
end
if ((~isempty(digestopt.varptm)) && (~isstruct(digestopt.varptm)))
        %||((~isempty(digestopt.fixedptm)) &&(~isstruct(digestopt.fixedptm))))
    error('MATLAB:GlycoPAT:INCORRECTTYPE','Wrong input type for the propertyname');
end
if dispprogress
    wtbar = waitbar(0, 'Performing protein digestion...');
    pause(0.0001);
end

% Reformat PTM information to be used in digestion
fixedptmopt = digestopt.fixedptm;
inputvarptmopt = digestopt.varptm;
numcores = digestopt.doparacomp;
varptmopt = [];
varptmuniind = 0;
if ~isempty(inputvarptmopt)
    tempvarptm = cell(length(inputvarptmopt),1);
    for ii = 1:length(tempvarptm)
        tempvarptm{ii} = [num2str(inputvarptmopt(ii).isglymod),...
            num2str(inputvarptmopt(ii).protpos),inputvarptmopt(ii).aaresidue];
    end
    [~,ind,varptmuniind] = unique(tempvarptm);
    varptmopt.aaresidue = {inputvarptmopt(ind).aaresidue};
    varptmopt.mod = {inputvarptmopt(ind).mod};
    varptmopt.protpos = {inputvarptmopt(ind).protpos};
    varptmopt.index = varptmuniind;
    varptmopt.isglymod = [inputvarptmopt(ind).isglymod];
    for ii = 1:length(varptmopt.mod)
        varptmopt.mod{ii} = ['{',num2str(ii),'}'];
    end
end

% Pack fixed PTM using variable PTM style - for easier processing
currentind = max(varptmuniind) + 1;
prevmod = {};
prevmodind = [];
if isstruct(fixedptmopt)
    for ii = 1:length(fixedptmopt)
        thisfixedptmnew = fixedptmopt(ii).new;
        [p,g,m] = breakGlyPep(thisfixedptmnew);
        inds = [];
        positions = [];
        if ~isempty(g)
            for jj = 1:length(g)
                thisstruct = g(jj).struct;
                structappeared = ismember(prevmod,thisstruct);
                if any(structappeared)
                    g(jj).struct = ['{',num2str(prevmodind(structappeared)),'}'];
                    inds = [inds,find(structappeared)];
                else
                    g(jj).struct = ['{',num2str(currentind),'}'];
                    prevmod = [prevmod,thisstruct];
                    prevmodind = [prevmodind,currentind];
                    inds = [inds,currentind];
                    currentind = currentind + 1;
                end
                positions = [positions,g(jj).pos];
            end
        end
        if ~isempty(m)
            for jj = 1:length(m)
                thisstruct = m(jj).struct;
                structappeared = ismember(prevmod,thisstruct);
                if any(structappeared)
                    m(jj).struct = ['{',num2str(prevmodind(structappeared)),'}'];
                    inds = [inds,find(structappeared)];
                else
                    m(jj).struct = ['{',num2str(currentind),'}'];
                    prevmod = [prevmod,thisstruct];
                    prevmodind = [prevmodind,currentind];
                    inds = [inds,currentind];
                    currentind = currentind + 1;
                end
                positions = [positions,m(jj).pos];
            end
        end
        [~,ind] = sort(positions);
        inds = inds(ind);
        fixedptmopt(ii).new_converted = joinGlyPep(p,g,m);
        fixedptmopt(ii).ind_converted = num2str(inds);
    end

    % Reserved - support for isobaric tags will be added in future release
    is_isotag = ismember({fixedptmopt.mod},Modification.isotagname);
    isobarictags = {fixedptmopt(is_isotag).mod};
    fixed.old = {fixedptmopt(~is_isotag).aaresidue};
    fixed.new = {fixedptmopt(~is_isotag).new_converted};
    fixed.protpos = {fixedptmopt(~is_isotag).protpos};
    
    % Here we use lower case letters to represent PTMs, this simplifies
    %     processing
    startChar = 97;  % placeholders start from char(97) = 'a'
    for ii = 1:length(fixed.old)
        fixed.replace{ii} = char(startChar);
        startChar = startChar + 1;
    end
    LastFixedChar = startChar - 1;  % This is the placeholder of the last fixed mod.
else
    startChar = 97;
    LastFixedChar = startChar - 1;
    fixed = 0;
    isobarictags = []; %fix later
end
varData = cell(1,4);
if(~isempty(varptmopt))
    varData{1} = varptmopt.mod;
    varData{2} = varptmopt.aaresidue;
    varData{3} = varptmopt.protpos;
    varData{4} = num2cell(varptmopt.isglymod);
else
    varData{1} = [];
    varData{2} = [];
    varData{3} = [];
    varData{4} = [];
end

% Process variable ptm
if isempty(varData{1})
    var.old={};
    var.letter={};
    var.number={};
    var.replaceWith={};
    modTable=[];
    modAacid=[];
    NXSTreplace='';
else
    % Consolidate varData into modAacid, modTable & var structures
    [var,modAacid,modTable,NXSTreplace]=...
        formatVariable(varData,startChar);
end

if dispprogress
    waitbar(0.05,wtbar);
end

% Now do the analysis, one protein at a time
if ischar(protseqheader)
    protseqheader = {protseqheader};
end
if numcores == 1
    ProdArray = cleaveproteins(protseq,protseqheader,fixed,enzname,...
        NXSTreplace,digestopt);
else
    parpool('local',numcores);
    protseq_dist = distributed(protseq);
    protseqheader_dist = distributed(protseqheader);
    spmd
        ProdArray_dist = cleaveproteins(getLocalPart(protseq_dist),...
            getLocalPart(protseqheader_dist),fixed,enzname,NXSTreplace,digestopt);
    end
    ProdArray = vertcat(ProdArray_dist{:});
end
nocleave = cellfun(@isempty,ProdArray);
ProdArray = ProdArray(~nocleave);

% Pause 0.1 sec to update waitbar
if dispprogress
    waitbar(0.1,wtbar);
end

if numcores == 1
    allProd = addvartopep(ProdArray,var,fixed,modAacid,modTable,NXSTreplace,LastFixedChar,digestopt);
    % Reserved - isobaric tags will be supported in future release
    if ~isempty(isobarictags)
        allProd = applyisobarictags(allProd,isobarictags);
    end
    if dispprogress
        waitbar(0.99,wtbar);
    end
else  % Parallel computing
    %% Method A: SPMD
    %     ProdArray_dist = distributed(ProdArray);
    %     spmd
    %         allProd_dist = addvartopep(getLocalPart(ProdArray_dist),var,fixed,modAacid,...
    %             modTable,NXSTreplace,LastFixedChar,digestopt);
    %         if ~isempty(isobarictags)
    %             allProd_dist = applyisobarictags(allProd_dist,isobarictags);
    %         end
    %     end
    %     allProd = vertcat(allProd_dist{:});
    %% Method A END
    
    %% Method B: asynchronous
    % 500 proteins per batch, these 500 proteins will be distributed to all
    %     cores
    prodchunksize = 500;
    numprod = ceil(length(ProdArray)/prodchunksize);
    ProdSto = cell(numprod,1);
    for ii = 1:numprod
        indstart = (ii-1) * prodchunksize + 1;
        indend = min(ii * prodchunksize,length(ProdArray));
        ProdSto{ii} = ProdArray(indstart:indend);
    end
    f(numprod) = parallel.FevalFuture;
    for ii = 1:numprod
        f(ii) = parfeval(@addvartopep,1,ProdSto{ii},var,fixed,modAacid,...
            modTable,NXSTreplace,LastFixedChar,digestopt);
    end
    updateWaitbarFuture = afterEach(f, @(~) waitbar(sum(strcmp('finished', {f.State}))/numel(f), wtbar), 1);
    allProd = fetchOutputs(f);
    afterAll(updateWaitbarFuture, @(wtbar) delete(wtbar), 0);
  %  close(wtbar);
    %% Method B END
end


finalallprod.glypep = allProd;
finalallprod.nocleave = protseqheader(nocleave);
finalallprod.FASTAhead = protseqheader(~nocleave);
finalallprod.enzyme = enzname;
finalallprod.missedmax = digestopt.missedmax;
finalallprod.minpeplen = digestopt.minpeplen;
finalallprod.maxpeplen = digestopt.maxpeplen;
finalallprod.minptm = digestopt.minptm;
finalallprod.maxptm = digestopt.maxptm;
finalallprod.varptm = digestopt.varptm;
finalallprod.protfilepath = digestopt.protfilepath;
for ii = 1:length(finalallprod.varptm)
    finalallprod.varptm(ii).varptmuniind = varptmuniind(ii);
end
finalallprod.fixedptm = fixedptmopt;
num_protein = num2str(length(protseq));  % Number of proteins digested
num_pep = num2str(sum(cellfun(@length,finalallprod.glypep)));  % Number of (glyco)peptide generated
reportstring = {['Protein digested: ' num_protein];...
    ['Generated peptide: ' num_pep]};
if ~isempty(finalallprod.nocleave)
    reportstring = [reportstring;'The following protein(s) was/were not digested:'];
    for ii = 1:length(finalallprod.nocleave)
        reportstring = [reportstring;finalallprod.nocleave{ii}];
    end
end
timeused = toc;
reportstring  = [reportstring;['Time used: ',num2str(timeused),' sec.']];
currenttime = clock;
reportstring  = [reportstring;['Finished at: ',num2str(currenttime(1)),'/',num2str(currenttime(2)),'/',num2str(currenttime(3)),'-'...
    ,num2str(currenttime(4)),':',num2str(currenttime(5)),':',num2str(floor(currenttime(6)))]];
outputstring = result2text(finalallprod,reportstring);  % Format results in pure text form
if numcores > 1
    delete(gcp('nocreate'));  % Close parallel computing environment
end
if dispprogress && numcores == 1
    close(wtbar);
end
end

function [ProdArray,counttot] = cleaveproteins(protseq,protseqheader,fixed,enzname,...
    NXSTreplace,digestopt)
% CLEAVEPROTEINS: perform cleavage on proteins
%
% Syntax:
% [ProdArray,counttot] = cleaveproteins(protseq,protseqheader,fixed,enzname,...
%     NXSTreplace,digestopt)
%
% Input:
% protseq: 1 x n cell array of strings. The amino acid sequences of
%     proteins.
% protseqheader: 1 x n cell array of strings. The FASTA header of
%     proteins.
% fixed: structure. Fixed PTM information. Fields are:
%     old: 1 x n cell array of strings. The amino acids to be modified.
%     new: 1 x n cell array of strings. The amino acids after modification.
%     protpos: 1 x n cell array of numbers. The position of these amino
%         acids.
%     replace: 1 x n cell array of strings. The single letter codes
%         representing modifications.
% enzname: 1 x n cell array of strings. The name of proteases.
% NXSTreplace: 1 x n cell array of strings. Special code for N-glycan
%     modification.
% digestopt: structure. Digestion options. See main function for detail.
%
% Output:
% ProdArray: n x 1 cell array of structures. Digested peptides. See
%     CLEAVEPROT for detail.
% counttot: Double. Number of generated peptides.
%
% Note:
% N/A
%
% Example:
% N/A
%
% Children function:
% N/A
%
% See Also:
% CLEAVEPROT

MinPepLen = digestopt.minpeplen;  % Min. length of digested peptide; default: 4
MaxPepLen = digestopt.maxpeplen;  % Max. length of digested peptide; default: 12
MissedMax = digestopt.missedmax;  % missed cleavage allowed; default: 0

nProt = length(protseqheader);
counttot=0;
ProdArray = cell(nProt,1);
FASTAhead = cell(nProt,1);

for np=1:nProt
    TempProt     = char(protseq(np));
    FASTAhead{np} = char(protseqheader(np));
    % Make NXST modifications since these are co-translated. They occur before fixed mod.
    % Sometimes there can be a site like 'NKT'. This would be N-linked but may get subsequently
    %     cleaved at the K site by Trypsin. Thus NST modification is done before proteolysis
    arra=regexp(TempProt,'N[ACDEFGHIKLMNQRSTVWY][ST]');
    if ~isempty(NXSTreplace)  % only done if N-glyans are present
        for kk = 1:length(arra)
            TempProt = [TempProt(1:arra(kk)-1),char(NXSTreplace),TempProt(arra(kk)+1:end)];
        end
    end
    if isstruct(fixed)
        for kk = 1:length(fixed.old)
            TempProt = regexprep(TempProt, char(fixed.old(kk)),...
                char(fixed.replace(kk))); % Apply fixed modifications
        end
    end
    ProdArray{np} = cleaveProt(char(TempProt),enzname,...
        MissedMax,MinPepLen,MaxPepLen); % calls the cleavage function
    if ~isempty(ProdArray{np})
        counttot = counttot + length(ProdArray{np}.pep);
    end
end
end

function allProd = addvartopep(ProdArray,var,fixed,modAacid,modTable,NXSTreplace,...
    LastFixedChar,digestopt)
% ADDVARTOPEP: combine variable and fixed PTM to transform digested peptide from
%     intermediate to final form.
%
% Syntax:
% allProd = addvartopep(ProdArray,var,fixed,modAacid,modTable,NXSTreplace,...
%     LastFixedChar,digestopt)
%
% Input:
% ProdArray: n x 1 cell array of structures. Digested peptides. See
%     CLEAVEPROT for detail.
% fixed: structure. Fixed PTMs. See CLEAVEPROTEINS for detail.
% var: structure. Variable PTMs. See FORMATVARIABLE for detail.
% modAacid: n x 1 cell array of strings. Identical to "var.old".
% modTable: n x 1 cell array of strings. Identical to "var.replaceWith".
% NXSTreplace: 1 x n cell array of strings. Special code for N-glycan
% LastFixedChar: double. Single letter code for the last fixed PTM.
% digestopt: structure. Digestion options. See DIGESTSGP for detail.
%
% Output:
% allProd: n x 1 cell array of m x 1 cell array of strings. Digested
%     glycopeptide in its final form.
%
% Note:
% N/A
%
% Example:
% N/A
%
% Children function:
% APPLYVARMOD
%
% See Also:
% DIGESTSGP  CLEAVEPROT  CLEAVEPROTEINS  FORMATVARIABLE  PERMUTATION

minPTM = digestopt.minptm;  % If the analysis is to be restricted to variable PTM only; default=0
maxPTM = digestopt.maxptm;  % Restrict number of variable PTM modifications; default=2
nProt = length(ProdArray);
allProd = cell(nProt,1);
for np=1:nProt
    Prod = ProdArray{np};
    % fixed modifications at site of digestion will inhibit cleavage
    % Variable modifications are applied on peptide products and not the original protein since this can
    %     restrict enzymatic digestion if the site of variable modification is
    %     also a site of protease digestion.
    % Here, the data available in 'var' is applied to modify all Prod
    finalProd = {};
    ProdOld = Prod.pep;
    Prod = applyVarMod(Prod,var);
    % Analysis is performed on each peptide Prod individually
    % Part A has to do with classifying all the modifications for a given
    %     peptide. This first part creates the Perm matrix that is used in Part B.
    % Here, it is ensured that no peptide has more than maxPTM modifications
    for i = 1:length(Prod.pep)
        originalPep = char(ProdOld(i));
        AA = cell2mat(Prod.pep(i));  % These two lines prevents variable mod. at the N-terminus (optional to include)
        str = [AA(1:length(AA)-1),originalPep(length(AA))];  % These two lines prevents variable mod. at the N-terminus (optional to include)
        aacidPos = find(double(str) > LastFixedChar);  % Find aacid positions with variable modifications
        aacidName = [];
        aacidNMod = [];
        for j = 1:length(aacidPos)
            aacidName = [aacidName,{str(aacidPos(j))}];  % edited by Gang LIu
            m = char(modAacid) == str(aacidPos(j));            % index for modification
            interim = find(~cellfun(@isempty,modTable(m,:)));
            aacidNMod = [aacidNMod, length(interim)];     % count how many modifications at each site
        end
        nPos = length(aacidPos);  % count how many modification sites
        Perm = [];
        if maxPTM == 0
            Perm = ones(nPos,1)';                 % return the original peptide
        elseif ((nPos <= maxPTM) && (nPos > 0))
            Perm = permutation(aacidNMod+1);      % +1 is added in order to include no modification condition
        elseif ((nPos > maxPTM) && (nPos > 0))
            c = combnk(1:nPos,maxPTM);            % Find all combinations of modifications
            sz1 = size(c);
            for k=1:sz1(1)
                nPosVec = zeros(1,maxPTM);
                for l=1:maxPTM
                    nPosVec(l)=aacidNMod(c(k,l));
                end
                YY=permutation(nPosVec+1);
                sz2=size(YY);
                XX = zeros(sz2(1),length(aacidNMod));
                for l=1:sz2(1)
                    XX(l,1:length(aacidNMod))=1;
                    for m=1:maxPTM
                        XX(l,c(k,m))=YY(l,m);
                    end
                end
                Perm = [Perm;XX];
            end
            Perm=unique(Perm,'rows');
        end
        % Part B uses data in Perm to replace amino acids in Prod.pep in order to
        % obtain the finalProd. Here it is ensured that no. of modifications >= minPTM
        if isempty(aacidPos)
            finalProd = [finalProd;str];
        else
            for j=1:size(Perm,1)
                nPTM=length(find(Perm(j,:)>1));
                if ((nPTM>=minPTM)||(all(Perm(j,:)==1)))
                    % make sure that the number of PTMs is more than minPTM, provided its not the original peptide
                    pepX='';
                    for k=1:length(str)
                        pos=find(aacidPos==k);
                        if(~isempty(pos))
                            modNumber=Perm(j,pos);
                            modAA =strcmp(str(k),modAacid);
                            if (modNumber==1)
                                intPTM='';
                            else
                                intPTM=modTable(modAA,modNumber-1);
                            end
                        else
                            intPTM='';
                        end
                        intProd = [originalPep(k),char(intPTM)];
                        pepX    = [pepX,intProd];
                    end
                    finalProd = [finalProd;pepX];
                end
            end
        end
    end
    % finally replace the fixed modifications and NST with the real ones
    if isstruct(fixed)
        for j=1:length(finalProd)
            for i=1:length(fixed.old)
                finalProd{j} = regexprep(char(finalProd(j)),fixed.replace(i),fixed.new(i));
            end
            finalProd{j}=(regexprep(char(finalProd(j)),NXSTreplace,'N'));
        end
    end
    finalProd        = unique(finalProd);
    allProd{np} = finalProd';
    %     if(dispprogress)
    %         waitbar(np/counttot,h);
    %     end
end
end


function [var,modAacid,modTable,NXSTreplace]=formatVariable(varData,startChar)
% FORMATVARIABLE: replace amino acids modified by variable PTM with single
%     letter code.
%
% Syntax:
% [var,modAacid,modTable,NXSTreplace]=formatVariable(varData,startChar)
%
% Input:
% vardata: 1 x 4 cell array. Elements are: 1. PTM strings; 2. amino acid
%     residues to be modified; 3. amino acid positions; 4. whether this PTM
%     is a glycan.
% startChar: double. The single letter used to describe the first PTM. This
%     value corresponds to ASCII if between 0 and 127, otherwise Unicode.
%
% Output:
% var: structure with 4 fields:
%     old: 1 x n cell array of strings. The single letter codes
%         representing PTMs.
%     replaceWith: 1 x n cell array of 1 x m cell array of strings. The
%         PTMs which the codes stand for.
%     letter: 1 x n cell array of strings. Target amino acid where PTMs
%         will be applied on.
%     number: 1 x n cell array of strings. Amino acid positions defined in
%         variable PTM input.
% modAacid: n x 1 cell array of strings. Identical to "var.old".
% modTable: n x 1 cell array of strings. Identical to "var.replaceWith".
% NXSTreplace: 1 x n cell array of strings. Special code for N-glycan
%     modification.
%
% Note:
% If a PTM has been given a specific modification site (defined in
%     varData{3}), output "var" will put this PTM in a separate slot
%
% Example:
% N/A
%
% Children function:
% N/A
%
% See Also:
% N/A

NXSTreplace = '';
var.old = {};
count = 1;

% First, deal with varData{2}, the amino acid characters
for i = 1:length(varData{2})
    aa = varData{2}{i};
    if ~(strcmp(aa,'X'))
        n = (size(aa,2)+1)/2;  % in case of multiple AAs can be mod. by 1 PTM, e.g. O-glycan
        for j = 1:n
            beenanalyzed = false;
            for k = 1:length(var.old)  % has this new PTM's AA been processed before?
                if strcmpi(aa(2*j-1),var.letter(k))
                    % if yes, for all AA associated, add it to the corresponding PTM list
                    var.replaceWith{k} = [var.replaceWith{k},varData{1}(i)];
                    beenanalyzed = true;
                end
            end
            if ~beenanalyzed  % new AA
                var.old(count) = {char(startChar)};
                var.replaceWith(count) = varData{1}(i);
                var.letter(count) = {aa(2*j-1)};
                var.number(count) = {0}; %%% edited by Gang
                startChar = startChar+1;
                firstChar = char(varData{1}(i));
                if ((char(var.letter(count)) == 'N') && (firstChar(1) == '{'))
                    NXSTreplace = var.old(count);  % special handling for N-glycan
                end
                count = count+1;  % Eventually "count" = number of unique AA to be modified
            end
        end
    end
end

% Second, deal with varData{3}, the position based modifications (amino acid
% based or position based)
for i = 1:length(varData{3})
    numaa  =  varData{3}{i};
    if( length(numaa)>1) ||( length(numaa) == 1 && numaa(1) ~= 0)
        for j = 1:length(numaa)
            beenanalyzed = false;
            for k = 1:length(var.old)
                varnum  =  var.number(k);
                if numaa(j) == varnum{:}
                    var.replaceWith{k} = [var.replaceWith{k},varData{1}(i)];
                    beenanalyzed  =  true;
                end
            end
            if (beenanalyzed  ==  false)
                var.old(count) = {(char(startChar))};
                var.replaceWith(count) = varData{1}(i);
                var.letter(count) = {'X'};
                var.number(count) = num2cell(numaa(j));
                startChar = startChar+1;
                count = count+1;
            end
        end
    end
end
modAacid = var.old';
for i = 1:length(var.replaceWith)
    if iscell(var.replaceWith{i})
        counter = length(var.replaceWith{i});
        for j = 1:counter
            modTable(i,j) = var.replaceWith{i}(j);
        end
    else
        modTable(i,1) = {var.replaceWith{i}};
    end
end
end

function Prod=applyVarMod(Prod,var)
% APPLYVARMOD: combine variable PTM to transform digested peptide from
%     intermediate to final form.
%
% Syntax:
% Prod=applyVarMod(Prod,var)
%
% Input:
% Prod: n x 1 cell array of structures. Digested peptides. See CLEAVEPROT
%     for detail.
% var: structure. Variable PTMs. See FORMATVARIABLE for detail.
%
% Output:
% Prod: n x 1 cell array of structures. Digested peptides with single
%     letter codes replaced by actual PTMs.
%
% Note:
% N/A
%
% Example:
% N/A
%
% Children function:
% N/A
%
% See Also:
% ADDVARTOPEP

for i=1:length(var.old)
    aa=char(var.letter(i));    % This allows the incorporation of variable modifications in comma-delimited form
    Numberaa=(length(aa)+1)/2;  % number of amino acids with identical modifications
    for j=1:Numberaa
        specificAA=aa(2*j-1);
        boo1=true;
        tempStr=char(var.replaceWith{i}(1));
        if ((specificAA=='N')&&(tempStr(1)=='{'))
            boo1=false;                                    % This is an N-glycan at N-X-S/T
        end
        if ((specificAA~='X')&&(boo1))    % if var.letter is not 'X' and this is not an NXS/T site
            varoldi = var.old(i);
            Prod.pep=regexprep(Prod.pep,specificAA,varoldi{:});  % modification at fixed aacid positions
        end
    end
    num=cell2mat(var.number(i));         % includes variable modifications when var.num ~=0
    if any(num~=0)  % Var. PTM with user defined modification site
        for j=1:length(num)
            for k=1:length(Prod.pep)
                if ((num(j)>=Prod.start(k))&&(num(j)<=Prod.fin(k)))
                    ProdTemp = Prod.pep(k);
                    ProdTemp = ProdTemp{:};
                    insertPos=num(j)-Prod.start(k);
                    ProdTemp=[ProdTemp(1:insertPos),cell2mat(var.old(i)),ProdTemp(insertPos+2:end)];
                    Prod.pep(k) = {ProdTemp};
                end
            end
        end
    end
end
end

function Perm=permutation(F)
% PERMUTATION: for an integer array of a given length, each member has
%     minimum of 1 and maximum given by the input, create all possible
%     combinations.
%
% Syntax:
% Perm=permutation(F)
%
% Input:
% F: 1 x n integer array. The maximum of each number.
%
% Output:
% Perm: m x n integer array. The combinations calculated from inputs. n
%     equals to the length of input "F".
%
% Note:
% N/A
%
% Example:
% Perm=permutation([1 2 3])
%
% Perm =
%
%      1     1     1
%      1     1     2
%      1     1     3
%      1     2     1
%      1     2     2
%      1     2     3
%
% Children function:
% N/A
%
% See Also:
% N/A

j = length(F);     % number of places
i = max(F);        % maximum value
Perm = varFor([],i,j);
a1 = size(Perm,1);
for k = 1:j
    l = 0;
    while  l < a1
        l = l+1;
        if Perm(l,k) > F(k)
            Perm(l,:) = [];
            l = l-1;
        end
        a1 = size(Perm,1);
    end
end
Perm = unique(Perm,'rows');
end

function outputstr = result2text(finalallprod,summarystring)
% RESULT2TEXT: rearrange the results into a format that can be written to a
%     text file directly.
%
% Syntax:
% outputstr = result2text(finalallprod,summarystring)
%
% Input:
% finalallprod: structure. Result of digestion. Fields are:
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%     Field name                        Format                                 Description
%     glypep                                n x 1 cell array of               Digested glycopeptides.
%                                                      m x 1 cell array of
%                                                      strings
%     nocleave                            1 x n cell array of               FASTA headers of proteins that generated no peptides.
%                                                      strings.
%     FASTAhead                       1 x n cell array of               FASTA headers of proteins in "glypep".
%                                                      strings.
%     enzyme                              1 x n cell array of               Name of the proteases.
%                                                      strings.
%     missedmax                        Double                                 Maximum number of missed cleavages allowed.
%     minpeplen                         Double                                 Mimimum length of digested peptides allowed.
%     maxpeplen                        Double                                 Maximum length of digested peptides allowed.
%     minptm                              Double                                 Minimum number of PTMs allowed on each peptide.
%     maxptm                             Double                                 Maximum number of PTMs allowed on each peptide.
%     varptm                               1 x n structure                    Variable PTMs. Fields are:  mod, aaresidue, protpos, isglymod.
%                                                                                                     See VARPTMREAD for detail.
%     fixedptm                            1 x n structure                    Fixed PTMs. Fields are:  aaresidue, new, protpos, mod.
%                                                                                                     See FIXEDPTMREAD for detail.
%     protfilepath                      String                                    The full path of protein .fasta file.
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%  summarystring: n x 1 cell array of strings. The summary of the
%      digestion, including the number of proteins digested, the number
%      of peptides generated and the time used.
%
% Output:
% outputstr: 1 x n string. Formatted string to be written to a text file.
%     Information has been catagorized and marked by headers starting with
%     "#".
%
% Note:
% Currently the headers are:
%     #SUMMARY
%     #VERSION
%     #PROTSEQPATH
%     #ENZYME
%     #PTM_FIXED
%     #PTM_VARIABLE
%     #PTM_MIN
%     #PTM_MAX
%     #PEPLEN_MIN
%     #PEPLEN_MAX
%     #MISSEDMAX
%     #PROTEIN_NOT_DIGESTED
%     #SEQUENCES
%     #END
% For the meaning of each header see DIGESTFILEANALY for detail.
%
% Example:
% N/A
%
% Children function:
% N/A
%
% See Also:
% DIGESTFILEANALY

outputstr = '';
%% SUMMARY
non_digest_txt_line = 3;  % Current version, only 2 lines.
outputstr = [outputstr,sprintf('%s\n','#SUMMARY')];
for ii = 1:non_digest_txt_line - 1
    tempstr = sprintf('%s\n',summarystring{ii});
    outputstr = [outputstr,tempstr];
end
%% TIME
outputstr = [outputstr,sprintf('%s\n',summarystring{end})];

outputstr = [outputstr,sprintf('%s\n','#VERSION')];
outputstr = [outputstr,sprintf('%s\n','1.1')];
outputstr = [outputstr,sprintf('%s\n','#PROTSEQPATH')];
outputstr = [outputstr,sprintf('%s\n',finalallprod.protfilepath)];

%% ENZYME
outputstr = [outputstr,sprintf('%s\n','#ENZYME')];
enzname = finalallprod.enzyme;
proteasedb = Protease.mklocaldb;
if(iscell(enzname))
    for k=1:length(enzname)
        if(proteasedb.isKey(upper(enzname{k})))
            ezncleaveexpr =  proteasedb(upper(enzname{k}));
            tempstr = sprintf('%i %s (name) %s (cleavage pattern)\n', k, ...
                enzname{k},ezncleaveexpr);
            outputstr = [outputstr,tempstr];
        else
            error('MATLAB:GlycoPAT:UNSUPPORTEDENZYME','UNSUPPORTED ENZYME TYPE');
        end
    end
else
    if(proteasedb.isKey(upper(enzname)))
        ezncleaveexpr =  proteasedb(upper(enzname));
        tempstr=sprintf('%s: (name) %s (cleavage pattern)\n', ...
            enzname,ezncleaveexpr);
        outputstr = [outputstr,tempstr];
    else
        error('MATLAB:GlycoPAT:UNSUPPORTEDENZYME','UNSUPPORTED ENZYME TYPE');
    end
end

%% FIXED PTM
outputstr = [outputstr,sprintf('%s\n','#PTM_FIXED')];
fixedptmopt = finalallprod.fixedptm;
if isstruct(fixedptmopt)
    for ii = 1 : length(fixedptmopt)
        tempstr=sprintf('%s %s %s %i \n',fixedptmopt(ii).ind_converted,...
            fixedptmopt(ii).aaresidue,fixedptmopt(ii).new,fixedptmopt(ii).protpos);
        outputstr = [outputstr,tempstr];
    end
end

%% VARIABLE PTM
outputstr = [outputstr,sprintf('%s\n','#PTM_VARIABLE')];
varptmopt = finalallprod.varptm;
if ~isempty(varptmopt)
    varptmindex = [finalallprod.varptm.varptmuniind];
    for ii = 1:length(varptmopt)
        tempstr=sprintf('%i %s %s %s\n',varptmindex(ii),varptmopt(ii).aaresidue,...
            varptmopt(ii).mod,num2str(varptmopt(ii).protpos));
        outputstr = [outputstr,tempstr];
    end
end

%% MINS & MAXS & MISSED
outputstr = [outputstr,sprintf('%s\n','#PTM_MIN'),sprintf('%i\n',finalallprod.minptm)];
outputstr = [outputstr,sprintf('%s\n','#PTM_MAX'),sprintf('%i\n',finalallprod.maxptm)];
outputstr = [outputstr,sprintf('%s\n','#PEPLEN_MIN'),sprintf('%i\n',finalallprod.minpeplen)];
outputstr = [outputstr,sprintf('%s\n','#PEPLEN_MAX'),sprintf('%i\n',finalallprod.maxpeplen)];
outputstr = [outputstr,sprintf('%s\n','#MISSEDMAX'),sprintf('%i\n',finalallprod.missedmax)];

%% NOT DIGESTED PROTEINS
outputstr = [outputstr,sprintf('%s\n','#PROTEIN_NOT_DIGESTED')];  % maybe add why they're not digested
if length(summarystring) - 1 > non_digest_txt_line
    for ii = non_digest_txt_line + 1 : length(summarystring) - 1
        tempstr=sprintf('%s\n',summarystring{ii});
        outputstr = [outputstr,tempstr];
    end
end

%% DIGESTED (GLY)PEP
outputstr = [outputstr,sprintf('%s\n','#SEQUENCES')];
glypeps = finalallprod.glypep;
FASTAhead = finalallprod.FASTAhead;
for ii =1:length(glypeps)
    tempstr    = sprintf('>%s\n',FASTAhead{ii});
    outputstr = [outputstr,tempstr];
    for j= 1:length(glypeps{ii})
        tempstr    = sprintf('%s\n',glypeps{ii}{j});
        outputstr = [outputstr,tempstr];
    end
end
%% FINISH
outputstr = [outputstr,sprintf('%s\n','#END')];
end