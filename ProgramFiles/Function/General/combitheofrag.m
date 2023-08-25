function finalfrag = combitheofrag(pepfrag,ptmfragsto,ptmfragstublen,pepseq,...
    ptmseq,ptmtype,ptmmass,gpcomp,fragnum,options)

% COMBITHEOFRAG: program for generating theoretical fragments through
% combination method
%
% Syntax:
% finalfrag = combitheofrag(pepfrag,ptmfragsto,ptmfragstublen,pepseq,...
% ptmseq,ptmtype,ptmmass,gpcomp,fragnum,options)
%
% Input:
% (Below "n" equals to number of all PTMs)
% pepfrag: structure, peptide fragments with "{...}" style markers representing PTMs.
% ptmfragsto: 1 x n cell array of structures. The storage place for PTM
%     fragments that corresponds to a certain fragmentation mode. Empty when
%     PTM is a non-glycan one.
% ptmfragstublen: 1 x n cell array of numbers. The length of glycan fragments.
% pepseq: string. Peptide sequence.
% ptmseq: 1 x n cell array of strings. Sequence of all PTMs.
% ptmtype: 1 x n numerical array. Type of each PTM. 1 means glycan, 2 means
%     non-glycan.
% ptmmass: 1 x n numerical array. Mass of each PTM.
% gpcomp: 1 x n numerical array. Glycopeptide composition. 1st number is
%     protein index, 2nd number is peptide index in this protein, 3rd and
%     onwards are indices of PTMs on this glycopeptide.
% fragnum: 1 x 3 numerical array. [npFrag, ngFrag, nmFrag]
% options: structure. Fields are:
%     fragmode: fragmentation mode.
%     stublen: maximum length of glycan residue on peptide fragments, only active
%         when fragmode is "HCD"
%     HCDstubfrag_override: if true, NgFrag for HCD glycan stub will be at
%         least 2. We recommend set this to TRUE.
%     HCDcofrag: if true, glycan residue size on peptide fragment will not
%         restrained by stublen. All possible fragments will be generated.
%     mode: 1 if you want to build sgp sequence for each fragment, 2 if you
%         don't. 2 is faster.
%     monosacislabile: when true, fragmentation of "f" and "s"  will not be
%         restricted by NgFrag.
%
% Output:
% finalfrag: structure, combined glycopeptide fragments.
%
% Note:
% N/A
%
% Example:
% load "SSAtestdata.mat" and run command:
% (Below "i", "j" can be differnet)
% theofrags = combitheofrag(combitf_input_pepfrag{i},ptmfragsto,ptmfragstublen,...
%     combitf_input_pepseq{i},ptmseq,ptmtype,ptmmass,combitf_input_gpcomp{i},...
%     combitf_input_fragnum{j},combitf_input_fragopt);
%
% Children function:
% FIXEDSUMDIST  COMBIPGMF

existingmonosac1let = 'hnsgfxzpukq';
maxstublen = options.stublen;
thisfragmode = options.fragmode;
switch options.mode
    case 1
        buildfragsgp = true;
    case 2
        buildfragsgp = false;
end
monosacislabile = options.monosacislabile;
simultaneousfrag = options.simultaneousfrag;
addoxoniumion = options.addoxoniumion;

[p,g,~] = breakGlyPep(pepseq);
ptmcompind = gpcomp(3:end);
for ii = 1:length(g)
    if gpcomp(2+ii) > 0
        g(ii).struct = ptmseq{gpcomp(2+ii)};
    else
        g(ii).struct = '';
    end
end
ptmcompind = ptmcompind(ptmcompind > 0);
finaloriginal = joinGlyPep(p,g,[]); % this is the sequence of the whole glycopeptide
[~,g,m] = breakGlyPep(finaloriginal);  % this is a simpler, less efficient approach
inpepvartype = ptmtype(ptmcompind);
if ~isempty(g)
    allgpos = [g.pos];  % this variable is for finding out which vptm is present on pepfrag
else
    allgpos = [];
end
if ~isempty(m)  % there is non-glycan vptm
    allmpos = [m.pos];
    modind = 1:length(m);
else  % there isn't
    allmpos = [];
end
ptmpos = sort([allgpos,allmpos]);
ptmcompind_available = ptmcompind(inpepvartype == 1 | inpepvartype == 2);

modunitindex = zeros(size(ptmcompind_available));  % this is a temporary variable containing ready to use unit index for non-glycan vptm
if ~isempty(allmpos)
    ind = 1;
    for ii = 1:length(inpepvartype)
        if inpepvartype(ii) == 2  % non-glycan vptm
            modunitindex(ptmcompind(ii)) = modind(ind);
            ind = ind + 1;
        end
    end
end

% NG/NM FRAG NUMBER DISTRIBUTION
ngfragmax = fragnum(2);
maxngfrags = cell(length(g),1);
mins = zeros(1,length(g));
maxs = ones(1,length(g)) * ngfragmax;
for ii = 1:length(maxngfrags)
    tempmaxngfrags = [];
    for jj = 0:ngfragmax
        tempmaxngfrags = [tempmaxngfrags;fixedsumdist(jj,mins(1:ii),maxs(1:ii))];
    end
    maxngfrags{ii} = tempmaxngfrags;
end
% MAXNGFRAGS: suppose 2 glycans on 1 peptide, ngfrag = 2, then glycan
% fragments can be [0,0], [0,1], [0,2], [1,0], [1,1], [2,0], they have sum <= 2.
% Retrieve these combinations by maxngfrags{k}, here k is the sum ngfrag
% numbers, e.g. k = 0 gets you [0,0], k = 1 gets you [1,0] and [0,1]

% nmfrag: on/not on peptide only, so it's either 0 or 1
nmfragmax = fragnum(3);
maxnmfrags = cell(length(m),1);
mins = zeros(1,length(m));
maxs = ones(1,length(m)) * nmfragmax;
for ii = 1:length(maxnmfrags)
    tempmaxnmfrags = [];
    for jj = 0:nmfragmax
        tempmaxnmfrags = [tempmaxnmfrags;fixedsumdist(jj,mins(1:ii),maxs(1:ii))];
    end
    badrow = sum(tempmaxnmfrags > 1,2) > 0;
    maxnmfrags{ii} = tempmaxnmfrags(~badrow,:);
end

ptmfragsto_f = ptmfragsto(ptmcompind_available);
ptmfragstublen_f = ptmfragstublen(ptmcompind_available);
ptmseq_f = ptmseq(ptmcompind_available);
ptmmass_f = ptmmass(ptmcompind_available);
ptmtype_f = ptmtype(ptmcompind_available);
% EXTRA LAYER OF SAFETY - JUST TO BE SAFE

%% PTM PRE-SORT BASED ON TYPE
frag_ptm_terminal = [];  % B-ions and middle ions (ngfrag >= 2, the part that's in the middle)
frag_ptm_intact = cell(size(ptmfragsto_f));  % 1 fragment only - glycan not fragmented.
frag_ptm_residue = cell(ngfragmax,length(ptmfragsto_f));  % Y-ions
fragstublen_ptm_residue = cell(ngfragmax,length(ptmfragsto_f));  % stub length of Y-ions
fragstublen_ptm_terminal = cell(ngfragmax,length(ptmfragsto_f));  % stub length of B-ions
fragstublen_ptm_middle = cell(ngfragmax,length(ptmfragsto_f));  % stub length of I-ions
glycanindoffset = 0;
%{
xyz = 0;
load 'C:\GlycoPAT\pfinal_test\11_01_2022 Syn4H Lib123 and Syn4H\Oglycoproteomics Syn4H\OSyn4H Bestisomers\Testvar.mat' xyz
if xyz >= 1225
    xyz = 0;
end
xyz = xyz + 1;
if xyz == 1225
    dbstop at 164 in combitheofrag
end
fprintf('%i \n',xyz)
save 'C:\GlycoPAT\pfinal_test\11_01_2022 Syn4H Lib123 and Syn4H\Oglycoproteomics Syn4H\OSyn4H Bestisomers\Testvar.mat' xyz
%}
for ii = 1:length(ptmfragsto_f)
    % i MARKS EACH PTM
    if ptmtype_f(ii) == 1  % if it's glycan
        tempfrag = ptmfragsto_f{ii};
        temptyp = {tempfrag.type};
        iontyp = zeros(size(tempfrag));
        for jj = 1:length(tempfrag)
            if ~isempty(strfind(temptyp{jj},'labile')) && strcmpi(temptyp{jj}(1:2),'-b')
                iontyp(jj) = 4;  % labile type b
            elseif ~isempty(strfind(temptyp{jj},'labile')) && strcmpi(temptyp{jj}(1:4),'none')
                iontyp(jj) = 5;  % labile intact
            elseif ~isempty(strfind(temptyp{jj},'labile')) && strcmpi(temptyp{jj}(1:2),'-y')
                iontyp(jj) = 6;  % labile type y
            elseif strcmpi(temptyp{jj}(1:2),'-b')
                iontyp(jj) = 0;
            elseif strcmpi(temptyp{jj},'none')
                iontyp(jj) = 1;
            elseif strcmpi(temptyp{jj}(1:2),'-y')
                iontyp(jj) = 2;
            elseif strcmpi(temptyp{jj},'ptm')
                iontyp(jj) = 3;
            end
            tempfrag(jj).unitindex{2} = tempfrag(jj).unitindex{2}(:,1) + glycanindoffset;
%             tempfrag(j).stublen = ptmfragstublen_f{i}(j);
        end
        glycanindoffset = glycanindoffset + size(tempfrag(iontyp == 1).unitindex{2},1);
        btype_ion = tempfrag(iontyp == 0);
        frag_ptm_intact{ii} = tempfrag(iontyp == 1 | iontyp == 3);
        frag_ptm_terminal = [frag_ptm_terminal,btype_ion];
        temp_ptm_residue = tempfrag(iontyp == 2);
        temp_ptmstublen_residue = ptmfragstublen_f{ii}(iontyp == 2);
        if monosacislabile
            frag_ptm_terminal = [frag_ptm_terminal,tempfrag(iontyp == 4)];
            temp_ptm_residue = [temp_ptm_residue,tempfrag(iontyp == 5 | iontyp == 6)];
            temp_ptmstublen_residue = [temp_ptmstublen_residue,ptmfragstublen_f{ii}(iontyp == 5 | iontyp == 6)];
        end
        temp_allngfragnum = [temp_ptm_residue.ngFrag];
        for jj = 1:max(temp_allngfragnum)
            frag_ptm_residue{jj,ii} = temp_ptm_residue(temp_allngfragnum == jj);
            fragstublen_ptm_residue{jj,ii} = temp_ptmstublen_residue(temp_allngfragnum == jj);
        end
        if ngfragmax > 0
            % SPECIAL TREATMENT: ADD AN EMPTY STRUCT TO SIMULATE Y0
            % NOT APPLICABLE WHEN NO GLYCAN FRAGMENTATION ALLOWED
            tempfrag = frag_ptm_intact{1,ii}(1);
            tempfrag.sgp = '';
            tempfrag.ngFrag = 1;
            tempfrag.mz = 18.0105647 + 1.007825032;
            tempfrag.type = ['-y-',tempfrag.original];
            tempfrag.unitindex{2} = zeros(0,1);
%             tempfrag.stublen = 0;
            frag_ptm_residue{1,ii} = [frag_ptm_residue{1,ii},tempfrag];
            fragstublen_ptm_residue{1,ii} = [fragstublen_ptm_residue{1,ii},0];
        end
    elseif ptmtype_f(ii) == 2  % if it's non-glycan PTM
        % AS OF NOW LEFT EMPTY
    end
end

%% PEPTIDE FRAGMENT SORTATION
inpepfragind_r = cell(size(pepfrag));  % WHICH PTM IS ON THIS PEP BACKBONE FRAGMENTS
for ii = 1:length(pepfrag)
    inpepfragind_r{ii} = find(ismember(ptmpos,pepfrag(ii).unitindex{1}));
end
woptmind = cellfun(@isempty,inpepfragind_r);
pep_wo_ptm = pepfrag(woptmind);  % PURE PEPTIDE FRAGMENT
[pep_wo_ptm.original] = deal(finaloriginal);
pep_w_ptm = pepfrag(~woptmind);  % PEPTIDE FRAG WITH PTM
pep_w_ptm_ind = inpepfragind_r(~woptmind);
[pep_w_ptm.original] = deal(finaloriginal);

if simultaneousfrag && fragnum(1) > 0 && fragnum(2) > 0  % GLY - PEP SIMULTANEOUS FRAGMENTATION
    % PUT LIMIT ON GLYCAN Y-FRAG, THEN COMBINE REGULARLY
    % A. all gly frag + intact pep
    combifragsA = combipgmf(pep_w_ptm(1),pep_w_ptm_ind(1),frag_ptm_intact,...
        frag_ptm_residue,ptmseq_f,ptmmass_f,ptmtype_f,...
        maxngfrags,maxnmfrags,modunitindex,buildfragsgp);
    if ~isempty(frag_ptm_residue) && maxstublen >= 0  % if maxstublen < 0, that's "unlimited stub length", no deletion necessary
        for ii = 1:size(frag_ptm_residue,1)
            % i MARKS EACH NGFRAG NUMBER
            for jj = 1:size(frag_ptm_residue,2)
                % j MARKS EACH VPTM
                if ~isempty(frag_ptm_residue{ii,jj})
                    goodind = fragstublen_ptm_residue{ii,jj} <= maxstublen;
                    frag_ptm_residue{ii,jj} = frag_ptm_residue{ii,jj}(goodind);
                end  % KEEP ONLY GLYCAN Y-IONS THAT'S SMALL ENOUGH
            end
        end
    end
    % stublen does not affect B-ions
    % B. truncated gly frag + pep frag
    combifragsB = combipgmf(pep_w_ptm(2:end),pep_w_ptm_ind(2:end),frag_ptm_intact,...
        frag_ptm_residue,ptmseq_f,ptmmass_f,ptmtype_f,...
        maxngfrags,maxnmfrags,modunitindex,buildfragsgp);
    combifrags = [combifragsA,combifragsB];
else  % no cofragmentation allowed
    % CAUTION: STUBLEN DOES NOT APPLY HERE
    % A. intact/labile + pepfrag
    tempmaxngfrags = maxngfrags;
    for ii = 1:length(tempmaxngfrags)
        templan = tempmaxngfrags{ii};
        tempmaxngfrags{ii} = templan(sum(templan,2) == 0,:);  % EXCLUDE NGFRAG > 0
    end
    combifragsA = combipgmf(pep_w_ptm(2:end),pep_w_ptm_ind(2:end),frag_ptm_intact,...
        frag_ptm_intact,ptmseq_f,ptmmass_f,ptmtype_f,...
        tempmaxngfrags,maxnmfrags,modunitindex,buildfragsgp);
    % (2:end), that means only pep frag. 1st element is combined with
    % glyfrag. This is to avoid duplicates.
    % Also use frag_ptm_intact as fake PTM fragments
    
    % B. glyfrag + intactpep
    if isempty(pep_w_ptm)
        pep_w_ptm = 0;
        pep_w_ptm_ind = 0;
    end
    combifragsB = combipgmf(pep_w_ptm(1),pep_w_ptm_ind(1),frag_ptm_intact,...
        frag_ptm_residue,ptmseq_f,ptmmass_f,ptmtype_f,...
        maxngfrags,maxnmfrags,modunitindex,buildfragsgp);
    % only limit peptide input, use full PTM frag reserves.
    combifrags = [combifragsB,combifragsA];
end

%% ASSEMBLE PTM - PEP
switch upper(thisfragmode)
    % MAKE NECESSARY TRIMMING ON GLYCAN & PEPTIDE FRAGMENTS
    case 'CID'
        % ONLY GLYCAN FRAGMENTATION:
        % B-ion + Y-ion (intact peptide backbone)
        % NOTHING NEEDS TO BE DONE HERE
    case 'HCD'
        % ADD OXONIUM ION AT THE END OF THE PROCESS
        if addoxoniumion && ~isempty(frag_ptm_terminal)
            frag_ptm_terminal = addoxoniumions(frag_ptm_terminal,existingmonosac1let);
        end
    case 'ETD'
        % NORMALLY NO GLYCAN FRAGMENTATION
        % c/z peptide ion
        % NOTHING NEEDS TO BE DONE HERE
    case 'ETCID'
        % BASICALLY IDENTICAL TO ETHCD
        % NO OXONIUM ION
        % NOTHING NEEDS TO BE DONE HERE
    case 'ETHCD'
        % B-ion + b/c/y/z peptide fragmentation allowed
        % ADD OXONIUM ION AT THE END OF THE PROCESS
        if addoxoniumion && ~isempty(frag_ptm_terminal)
            frag_ptm_terminal = addoxoniumions(frag_ptm_terminal,existingmonosac1let);
        end
    otherwise
        if addoxoniumion && ~isempty(frag_ptm_terminal)
            frag_ptm_terminal = addoxoniumions(frag_ptm_terminal,existingmonosac1let);
        end
end
finalfrag = [];
if ~isempty(frag_ptm_terminal)
    [frag_ptm_terminal.original] = deal(finaloriginal);
    finalfrag = [finalfrag,frag_ptm_terminal];
end
if ~isempty(pep_wo_ptm)
    finalfrag = [finalfrag,pep_wo_ptm];
end
if ~isempty(combifrags)
    finalfrag = [finalfrag,combifrags];
end
% Bring "none" to the first place - makes later steps easier
noneion = find(strcmpi({finalfrag.type},'none'));
noneion = noneion(1);
newseq = [noneion,setdiff(1:length(finalfrag),noneion)];
finalfrag = finalfrag(newseq);
end

function frag_vptm_terminal = addoxoniumions(frag_vptm_terminal,existingmonosac1let)
% ADDOXONIUMIONS: for HCD and ETHCD, add oxonium ions
%
% Syntax:
% frag_vptm_terminal = addoxoniumions(frag_vptm_terminal,existingmonosac1let)
%
% Input:
% frag_vptm_terminal: B-type fragment ions.
% existingmonosac1let: string, 1 letter monosac. available now.
% Output:
% frag_vptm_terminal: B-type fragment ions with additional oxonium ions.
%
% Note:
% This function searches B-type ions and see if "s", "h" and "n" is present
%     at terminal of structure. If so, the corresponding oxonium ions will be
%     added.
%
% Example:
% N/A
%
% Children function:
% N/A

terminalfragsgps = {frag_vptm_terminal.sgp};

borrowunitind = find(strcmpi(terminalfragsgps,'{s}'));
buildparent = false;
if isempty(borrowunitind)
    borrowunitind = find(~cellfun(@isempty,strfind(terminalfragsgps,'{s')),1,'first');
    buildparent = true;
end
if any(borrowunitind)
    borrowunitind = borrowunitind(1);
    len = length(terminalfragsgps);
    frag_vptm_terminal(len+1).original = frag_vptm_terminal(borrowunitind).original;
    frag_vptm_terminal(len+1).sgp = '{s-water}';
    frag_vptm_terminal(len+1).nmFrag = 0;
    frag_vptm_terminal(len+1).npFrag = 0;
    frag_vptm_terminal(len+1).ngFrag = 1;
    frag_vptm_terminal(len+1).mz = 274.092675035;
    frag_vptm_terminal(len+1).type = '-b{s-water}oxo';
    frag_vptm_terminal(len+1).charge = 1;
%     frag_vptm_terminal(len+1).stublen = 1;
    tempunitindex = cell(1,3);
    borrowedsgp = frag_vptm_terminal(borrowunitind).sgp;
    borrowedunitindex2 = frag_vptm_terminal(borrowunitind).unitindex{2};
    tgtsgppos = strfind(borrowedsgp,'s');
    monosacpos = regexp(borrowedsgp,['[',existingmonosac1let,']']);
    tempunitindex{2} = borrowedunitindex2(monosacpos == tgtsgppos(end));
    frag_vptm_terminal(len+1).unitindex = tempunitindex;
    if buildparent
        frag_vptm_terminal(len+2).original = frag_vptm_terminal(borrowunitind).original;
        frag_vptm_terminal(len+2).sgp = '{s}';
        frag_vptm_terminal(len+2).nmFrag = 0;
        frag_vptm_terminal(len+2).npFrag = 0;
        frag_vptm_terminal(len+2).ngFrag = 2;
        frag_vptm_terminal(len+2).mz = 292.103240032;
        frag_vptm_terminal(len+2).type = '-b{s}';
        frag_vptm_terminal(len+2).charge = 1;
        frag_vptm_terminal(len+2).unitindex = tempunitindex;
%         frag_vptm_terminal(len+2).stublen = 1;
    end
end

borrowunitind = find(strcmpi(terminalfragsgps,'{h}'));
buildparent = false;
if isempty(borrowunitind)
    borrowunitind = find(~cellfun(@isempty,strfind(terminalfragsgps,'{h')));
    buildparent = true;
end
if any(borrowunitind)
    borrowunitind = borrowunitind(1);
    len = length(frag_vptm_terminal);
    frag_vptm_terminal(len+1).original = frag_vptm_terminal(borrowunitind).original;
    frag_vptm_terminal(len+1).sgp = '{h-water}';
    frag_vptm_terminal(len+1).nmFrag = 0;
    frag_vptm_terminal(len+1).npFrag = 0;
    frag_vptm_terminal(len+1).ngFrag = 1;
    frag_vptm_terminal(len+1).mz = 145.050082435;
    frag_vptm_terminal(len+1).type = '-b{h-water}oxo';
    frag_vptm_terminal(len+1).charge = 1;
%     frag_vptm_terminal(len+1).stublen = 1;
    tempunitindex = cell(1,3);
    borrowedsgp = frag_vptm_terminal(borrowunitind).sgp;
    borrowedunitindex2 = frag_vptm_terminal(borrowunitind).unitindex{2};
    tgtsgppos = strfind(borrowedsgp,'h');
    monosacpos = regexp(borrowedsgp,['[',existingmonosac1let,']']);
    tempunitindex{2} = borrowedunitindex2(monosacpos == tgtsgppos(end));
    frag_vptm_terminal(len+1).unitindex = tempunitindex;
    if buildparent
        frag_vptm_terminal(len+2).original = frag_vptm_terminal(borrowunitind).original;
        frag_vptm_terminal(len+2).sgp = '{h}';
        frag_vptm_terminal(len+2).nmFrag = 0;
        frag_vptm_terminal(len+2).npFrag = 0;
        frag_vptm_terminal(len+2).ngFrag = 2;
        frag_vptm_terminal(len+2).mz = 163.060647132;
        frag_vptm_terminal(len+2).type = '-b{h}';
        frag_vptm_terminal(len+2).charge = 1;
        frag_vptm_terminal(len+2).unitindex = tempunitindex;
%         frag_vptm_terminal(len+2).stublen = 1;
    end
end

borrowunitind = find(strcmpi(terminalfragsgps,'{n}'));
buildparent = false;
if isempty(borrowunitind)
    borrowunitind = find(~cellfun(@isempty,strfind(terminalfragsgps,'{n')));
    buildparent = true;
end
if any(borrowunitind)
    borrowunitind = borrowunitind(1);
    len = length(frag_vptm_terminal);
    frag_vptm_terminal(len+1).original = frag_vptm_terminal(borrowunitind).original;
    frag_vptm_terminal(len+1).sgp = '{n-water}';
    frag_vptm_terminal(len+1).nmFrag = 0;
    frag_vptm_terminal(len+1).npFrag = 0;
    frag_vptm_terminal(len+1).ngFrag = 1;
    frag_vptm_terminal(len+1).mz = 186.076631035;
    frag_vptm_terminal(len+1).type = '-b{n-water}oxo';
    frag_vptm_terminal(len+1).charge = 1;
%     frag_vptm_terminal(len+1).stublen = 1;
    tempunitindex = cell(1,3);
    borrowedsgp = frag_vptm_terminal(borrowunitind).sgp;
    borrowedunitindex2 = frag_vptm_terminal(borrowunitind).unitindex{2};
    tgtsgppos = strfind(borrowedsgp,'n');
    monosacpos = regexp(borrowedsgp,['[',existingmonosac1let,']']);
    tempunitindex{2} = borrowedunitindex2(monosacpos == tgtsgppos(end));
    frag_vptm_terminal(len+1).unitindex = tempunitindex;
    frag_vptm_terminal(len+2).original = frag_vptm_terminal(borrowunitind).original;
    frag_vptm_terminal(len+2).sgp = '{n}';
    frag_vptm_terminal(len+2).nmFrag = 0;
    frag_vptm_terminal(len+2).npFrag = 0;
    frag_vptm_terminal(len+2).ngFrag = 1;
    frag_vptm_terminal(len+2).mz = 168.066066035;
    frag_vptm_terminal(len+2).type = '-b{n-water-water}oxo';
    frag_vptm_terminal(len+2).charge = 1;
    frag_vptm_terminal(len+2).unitindex = tempunitindex;
%     frag_vptm_terminal(len+2).stublen = 1;
    frag_vptm_terminal(len+3).original = frag_vptm_terminal(borrowunitind).original;
    frag_vptm_terminal(len+3).sgp = '{n}';
    frag_vptm_terminal(len+3).nmFrag = 0;
    frag_vptm_terminal(len+3).npFrag = 0;
    frag_vptm_terminal(len+3).ngFrag = 1;
    frag_vptm_terminal(len+3).mz = 144.066066335;
    frag_vptm_terminal(len+3).type = '-b{n-C2H4O2}oxo';
    frag_vptm_terminal(len+3).charge = 1;
    frag_vptm_terminal(len+3).unitindex = tempunitindex;
%     frag_vptm_terminal(len+3).stublen = 1;
    frag_vptm_terminal(len+4).original = frag_vptm_terminal(borrowunitind).original;
    frag_vptm_terminal(len+4).sgp = '{n}';
    frag_vptm_terminal(len+4).nmFrag = 0;
    frag_vptm_terminal(len+4).npFrag = 0;
    frag_vptm_terminal(len+4).ngFrag = 1;
    frag_vptm_terminal(len+4).mz = 138.055501335;
    frag_vptm_terminal(len+4).type = '-b{n-CH6O3}oxo';
    frag_vptm_terminal(len+4).charge = 1;
    frag_vptm_terminal(len+4).unitindex = tempunitindex;
%     frag_vptm_terminal(len+4).stublen = 1;
    if buildparent
        frag_vptm_terminal(len+5).original = frag_vptm_terminal(borrowunitind).original;
        frag_vptm_terminal(len+5).sgp = '{n}';
        frag_vptm_terminal(len+5).nmFrag = 0;
        frag_vptm_terminal(len+5).npFrag = 0;
        frag_vptm_terminal(len+5).ngFrag = 2;
        frag_vptm_terminal(len+5).mz = 204.087196032;
        frag_vptm_terminal(len+5).type = '-b{n}';
        frag_vptm_terminal(len+5).charge = 1;
        frag_vptm_terminal(len+5).unitindex = tempunitindex;
%         frag_vptm_terminal(len+5).stublen = 1;
    end
end
end