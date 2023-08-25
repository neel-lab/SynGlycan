function [mass,comp] = assemblegp(prot,ptmmass,ptmmassfix,ptminds,ptmstruct,options)
% ASSEMBLEGP: combine peptide backbone and PTMs, return mass and
%     composition of the product
% 
% Syntax:
% [mass,comp] = assemblegp(prot,ptmmass,ptmmassfix,ptminds,ptmstruct,options)
% 
% Input:
% prot: n x 1 cell array of m x 1 cell array of strings. The peptide 
%     backbones of candidate glycopeptides. n equals to the number
%     of proteins, m equals to the number of glycopeptide's backbone.
%     The  sequence of these backbones should be in the form of
%     "PEPT{i}DE", in which "PEPTDE" are the amino acids, "{i}" is
%     the locator of PTMs on this peptide, "i" is an integer. Only PTMs 
%     with an identical loca tor may be placed on this site.
% ptmmass: n x 1 double. The mass of each PTM.
% ptmmassfix: n x 1 double. The mass fix of each PTM to compensate mass
%     changes when binding to the peptide backbone, for example water loss.
% ptminds: n x 1 double. The locator of each PTM. Each PTM has one locator.
% ptmstruct: n x 1 cell array of strings. PTM sequence.
% options: currently is available field is:
%     "mode": 1 or 2. In mode 1, a locator must be occupied by a
%     corresponding PTM. In mode 2, this site can be left vacant.
%     Mode 2 is reserved for a possible future improvement which aims to
%     reduce digestion file size.
% 
% Output:
% mass: the theoretical monosiotopic mass of all combined glycopeptides.
% comp: n x 1 cell array, the composition of all combined glycopeptides.
%     1st number is protein number, 2nd is peptide number in this protein,
%     3rd and onwards are PTM numbers. Here the PTM numbers are not
%     locators, they denote the location in input "ptmmass", "ptmmassfix",
%     "ptminds" and "ptmstruct". See also FINDMS1MATCH.
% 
% Note:
% As of now the use of mode 2 is not recommended as it will create
%     duplicated glycopeptide combinations.
%
% Example:
% N/A
%
% Children function: 
% N/A

%% NUMBER OF VAR PTM
[~,~,ind] = unique(ptminds);
ptmspecies = cell(max(ind),1);
for ii = 1:max(ind)
    tempind = ind == ii;
    ptmspecies{ii} = ptmstruct(tempind);
end
numptmspecies = cellfun(@length,ptmspecies);  % how many mods each index number has
ptmrealidfix = zeros(1,length(numptmspecies));
% this fix value will help retrieve the varptm frag from storage
% each varptm type has members, use this fix value we can have:
% species x -> fix value for specie x + yth member = which line to retrieve
% from storage
if length(numptmspecies) > 1
    for ii = 2:length(numptmspecies)
        ptmrealidfix(ii) = ptmrealidfix(ii-1) + numptmspecies(ii-1);
    end
end

%% try to save some memory - abandoned because may not be that efficient
% pepnum = cellfun(@length,prot);
% writeind = 1;
% allpep = cell(sum(pepnum),1);
% allpepind = zeros(sum(pepnum),2);
% for i = 1:length(prot)
%     thisprot = prot{i};
%     allpep(writeind:writeind + length(thisprot) - 1) = thisprot(:);
%     allpepind(writeind:writeind + length(thisprot) - 1,1) = ones(length(thisprot),1)*i;
%     allpepind(writeind:writeind + length(thisprot) - 1,2) = (1:length(thisprot))';
%     writeind = writeind + length(thisprot);
% end
% [~,~,tempunipepind] = unique(allpep,'stable');
% unipepind = cell(size(prot));
% for i = 1:length(prot
% clear allpep;  % release memory

%% ASSEMBLE SGP
% GlycoPAT 1 design: if a glycopeptide can be modified by varmod in
% different ways (e.g. only 1 N-glycan, but attached to different
% Asn in pep), there will be more than 1 pep, showing the
% difference in N-glycosylation
gpcounter = 0;
for ii = 1:length(prot)
    thisprot = prot{ii};
    for jj = 1:length(thisprot)
        thispep = thisprot{jj};
        thisptmstruct = regexp(thispep,'{[0-9]+}','match');
        switch options.mode
            case 1  % GLYCOPAT 1 DIGEST DESIGN: IF A PTM SITE IS SHOWN OCCUPIED
                %  THEN THIS SITE MUST HAVE SOMETHING ON IT
                if ~isempty(thisptmstruct)
                    ptmspeciesnum = zeros(size(thisptmstruct));
                    for kk = 1:length(thisptmstruct)
                        ptmind = str2double(thisptmstruct{kk}(2:end-1));
                        ptmspeciesnum(kk) = numptmspecies(ptmind);
                    end
                    gpcounter = gpcounter + prod(ptmspeciesnum);
                else
                    gpcounter = gpcounter + 1;
                end
            case 2
                % DON'T USE NOW! IN MODE 2 MARKED PTM SITE MAY BE
                % VACANT
                if ~isempty(thisptmstruct)
                    ptmspeciesnum = zeros(size(thisptmstruct));
                    for kk = 1:length(thisptmstruct)
                        ptmind = str2double(thisptmstruct{kk}(2:end-1));
                        ptmspeciesnum(kk) = numptmspecies(ptmind) + 1;  % +1 for no-PTM
                    end
                    gpcounter = gpcounter + prod(ptmspeciesnum);
                else
                    gpcounter = gpcounter + 1;
                end
        end
    end
end
mass = zeros(gpcounter,1);
comp = cell(gpcounter,1);
writeind = 1;
for ii = 1:length(prot)
    thisprot = prot{ii};
    for jj = 1:length(thisprot)
        thispep = thisprot{jj};
        [thisptmpos,thisptmstruct] = regexp(thispep,'{[0-9]+}','start','match');
        if ~isempty(thisptmpos)
            temppep = regexprep(thispep,'{[0-9]+}','');
            pepmass = pepMW(temppep,1);
            whichptm = zeros(size(thisptmpos));
            howmanyptm = zeros(size(thisptmpos));
            % for each varptm represented by a number how many
            % members are there?
            % each varptm is independent from each other
            for kk = 1:length(howmanyptm)
                whichptm(kk) = str2double(thisptmstruct{kk}(2:end-1));
                howmanyptm(kk) = numptmspecies(whichptm(kk));
            end
            thisptmrealidfix = ptmrealidfix(whichptm);
            switch options.mode
                case 1
                    thisptmmodcombi = getallcombi(howmanyptm);
                case 2
                    % DON'T USE NOW! SHOULD BE MORE COMPLICATED THAN THIS
                    % CURRENT DESIGN CAN'T DISTINGUISH FIXED FROM VARIABLE
                    % PTM
                    thisptmmodcombi = getallcombi(howmanyptm + 1) - 1;
            end
            thisgpmasses = zeros(size(thisptmmodcombi,1),1);
            thiscompes = cell(size(thisptmmodcombi,1),1);
            for kk = 1:size(thisptmmodcombi,1)
                tempcompes = thisptmmodcombi(kk,:) + thisptmrealidfix;
                % these will be the real ptm index numbers, use
                % directly to retrieve pre-fragmented ptm
                thisgpmasses(kk) = pepmass;
                for ll = 1:length(tempcompes)
                    if tempcompes(ll) > 0
                        thisgpmasses(kk) = thisgpmasses(kk) + ptmmass(tempcompes(ll)) + ptmmassfix(tempcompes(ll));
                    end
                end
                thiscompes{kk} = [ii,jj,tempcompes];
            end
            mass(writeind:writeind + kk - 1) = thisgpmasses;
            comp(writeind:writeind + kk - 1) = thiscompes;
            writeind = writeind + kk;
        else
            mass(writeind) = pepMW(thispep,1);
            comp{writeind} = [ii,jj];
            writeind = writeind + 1;
        end
    end
end
end