function combifrags = combipgmf(pepfrag,inpepfragind_r,frag_ptm_intact,...
    frag_ptm_residue,ptmseq_f,ptmmass_f,ptmtype_f,ngfragdist,nmfragdist,...
    modunitindex,buildfragsgp)
% COMBIPGMF: generate theoretical fragments through combination (those with peptides)
%
% Syntax:
% combifrags = combipgmf(pepfrag,inpepfragind_r,frag_ptm_intact,...
% frag_ptm_residue,ptmseq_f,ptmmass_f,ptmtype_f,ngfragdist,nmfragdist,...
% modunitindex,buildfragsgp)
%
% Input:
% (Below "n" equals to the number of PTMs that appears on this glycopeptide)
% pepfrag: structure, peptide fragments. Each one must have PTM attached.
% inpepfragind_r: 1 x n cell array of numerical arrays. Indices of PTMs
%     attached to each peptide fragment. ("_r" means "refined")
% frag_ptm_intact: 1 x n cell array of structures. The intact form of PTM
%     fragments.
% frag_ptm_residue: 1 x n cell array of structures. The fragmented form of
%     PTM fragments.
% ptmseq_f: 1 x n cell array of strings. The sequence of PTMs on this
%     peptide. ("_f" means frag mode specific)
% ptmmass_f: 1 x n cell array of strings. The mass of PTMs on this
%     peptide.
% ptmtype_f: 1 x n cell array of strings. The type of PTMs on this
%     peptide. 1 for glycan, 2 for non-glycan.
% ngfragdist: 1 x m cell array. The glycan cleavage number distribution
%     plans. The Mth cell contains the plan to use when there are M
%     glycans. "Plan" means the NgFrag to be applied on specific glycans.
%     For each plan, the sum of all NgFrag numbers equal to the max NgFrag
%     specified when starting the parent program.
% nmfragdist: 1 x m cell array. The non-glycan PTM cleavage number
%     distribution plans.
% modunitindex: 1 x n numerical array. The unit indices for non-glycan PTMs.
% buildfragsgp: true/false, build SGP sequence for fragment or not.
% 
% Output:
% combifrags: 1 x n structure. The fragments of this glycopeptide
%
% Note:
% N/A
%
% Example:
% load "testdata_combipgmf.mat" and run command:
% combifrags = combipgmf(pepfrag,inpepfragind_r,frag_ptm_intact,...
%     frag_ptm_residue,ptmseq_f,ptmmass_f,ptmtype_f,ngfragdist,nmfragdist,...
%     modunitindex,buildfragsgp)
%
% Children function:
% N/A

%% ESTIMATE MAXIMUM NUMBER OF COMBINED FRAGMENTS - 
% PREALLOCATE MEMORY SPACE FOR BETTER SPEED
estifraglen = 0;
ptmfragnum = sum(cellfun(@length,frag_ptm_residue),1);
for i = 1:length(inpepfragind_r)
    estifraglen = estifraglen + prod(ptmfragnum(inpepfragind_r{i}) + 1 +...
        (ptmtype_f(inpepfragind_r{i}) == 2));
end

% pre-fetch unit indices (AA & PTM) because they are used heavily,
% this saves time rather than retrieving from data structure everytime.
frag_ptm_intact_glycanind = cell(size(frag_ptm_intact));
for i = 1:length(frag_ptm_intact)
    if ~isempty(frag_ptm_intact{i})
        frag_ptm_intact_glycanind{i}{1} = frag_ptm_intact{i}.unitindex{2}(:,1);
    end
end

frag_ptm_residue_glycanind = cell(size(frag_ptm_residue));
for i = 1:size(frag_ptm_residue,1)
    for j = 1:size(frag_ptm_residue,2)
        tempfrag_ptm_residue = frag_ptm_residue{i,j};
        tempfrag_ptm_residue_glycanind = cell(size(tempfrag_ptm_residue));
        for k = 1:length(tempfrag_ptm_residue)
            tempfrag_ptm_residue_glycanind{k} = tempfrag_ptm_residue(k).unitindex{2}(:,1);
        end
        frag_ptm_residue_glycanind{i,j} = tempfrag_ptm_residue_glycanind;
    end
end

original = cell(estifraglen,1);
sgp = cell(estifraglen,1);
npFrag = cell(estifraglen,1);
ngFrag = cell(estifraglen,1);
nmFrag = cell(estifraglen,1);
mz = cell(estifraglen,1);
type = cell(estifraglen,1);
charge = cell(estifraglen,1);
unitindex = cell(estifraglen,1);
stublen = cell(estifraglen,1);
writeind = 1;
for i = 1:length(pepfrag)
    % i MARKS PEPTIDE BACKBONE
    thisinpepfragind_r = inpepfragind_r{i};  % which PTM are on this pep frag
    thisptmtyp = ptmtype_f(thisinpepfragind_r);  % these numbers are in the order of their appearance on pep already
    numgly = sum(thisptmtyp == 1);
    numnongly = sum(thisptmtyp == 2);
    basestruct = pepfrag(i);  % base, or peptide backbone, is the foundation for combining gp frag
    baseoriginal = basestruct.original;
    basesgp = basestruct.sgp;
    basemz = basestruct.mz;
    baseunitindex = basestruct.unitindex;
    basenpFrag = basestruct.npFrag;
    basetype = basestruct.type;
    basestublen = [];
    if numnongly > 0  % there is non-glycan PTM
        if numgly > 0  % glycan and non-glycan
            this_nxfrag_plan = combigroups(ngfragdist{numgly},...
                nmfragdist{numnongly});  % distribution of fragmentation numbers
            ngnmfragnum = combigroups(sum(ngfragdist{numgly},2),...
                sum(nmfragdist{numnongly},2));  % total number of fragments, to be written into final product
        else  % non-glycan only
            this_nxfrag_plan = nmfragdist{numnongly};
            ngnmfragnum = zeros(size(this_nxfrag_plan,1),2);
            ngnmfragnum(:,2) = sum(this_nxfrag_plan,2);
        end
    else  % only glycan
        this_nxfrag_plan = ngfragdist{numgly};
        ngnmfragnum = zeros(size(this_nxfrag_plan,1),2);
        ngnmfragnum(:,1) = sum(this_nxfrag_plan,2);
    end
    [~,ind] = sort(thisptmtyp);  % originally it's [ng1,ng2,...,nm1,nm2,...]
    [~,ind] = sort(ind);  % ind must be resorted before using, now ind at each position means 
    % which element to be placed here
    this_nxfrag_plan = this_nxfrag_plan(:,ind);  % using this "ind", the plan has ng and nmfrag
    % distributed according to their relative position on peptide
    % each row is a nm/ngfrag distribution pattern. Use these numbers to
    % find appropriate PTM fragment to plant on pep
    
    % HERE "PLAN" IS FRAG NUMBER DISTRIBUTION
    if buildfragsgp
        [p,g,m] = breakGlyPep(basesgp);
    end
    for j = 1:size(this_nxfrag_plan,1)
        % j MARKS NUMFRAG DISTRIBUTION
        tempptm_frag = cell(size(thisinpepfragind_r));
        tempptm_frag_glycanind = cell(size(thisinpepfragind_r));
        tempptm_typ = ptmtype_f(thisinpepfragind_r);
        tempptm_seq = ptmseq_f(thisinpepfragind_r);
        tempptm_mass = ptmmass_f(thisinpepfragind_r);
        tempptm_stublen = {};
        % assemble to-be-used ptm fragments
        for k = 1:length(tempptm_frag)
            % k MARKS EACH PTM
            if this_nxfrag_plan(j,k) == 0  % use intact form of glycan OR a non-glycan PTM IS here
                % Here 0 means no fragmentation happened
                if tempptm_typ(k) == 1  % this is an intact glycan
                    tempptm_frag{k} = frag_ptm_intact{thisinpepfragind_r(k)};
                    tempptm_frag_glycanind{k} = frag_ptm_intact_glycanind{thisinpepfragind_r(k)};
                    tempptm_stublen = [tempptm_stublen,frag_ptm_intact{thisinpepfragind_r(k)}.stublen];
                elseif tempptm_typ(k) == 2  % this is a non-glycan
                    tempptm_frag{k} = 1;
                end
            else  % use fragmented form of glycan OR non-glycan PTM is NOT here
                if tempptm_typ(k) == 1 && size(frag_ptm_residue,1) >= this_nxfrag_plan(j,k)
                    % the 2nd part puts a limitation on ngFrag number: when
                    % available glycan frag supports only ngFrag = 1, if
                    % plan has ngFrag = 2, nothing will happen. This leaves
                    % an empty space in the variable, which will be
                    % rejected, no SGP building will take place
                    tempptm_frag{k} = frag_ptm_residue{this_nxfrag_plan(j,k),thisinpepfragind_r(k)};
                    tempptm_frag_glycanind{k} = frag_ptm_residue_glycanind{this_nxfrag_plan(j,k),...
                        thisinpepfragind_r(k)};
                    tempptm_stublen = [tempptm_stublen,...
                        [frag_ptm_residue{this_nxfrag_plan(j,k),thisinpepfragind_r(k)}.stublen]];
                elseif tempptm_typ(k) == 2  % this is a non-glycan
                    tempptm_frag{k} = 0;  % 0 occupies the space so it's not empty, so fragment will be built
                end
            end
        end
        asmbplan = getallcombi(cellfun(@length,tempptm_frag));
        
        % "asmbplan" is the receipe for fragment assembly
        % ASSEMBLE START
        for k = 1:size(asmbplan,1)
            % k MARKS EACH ASSEMBLY PLAN (which fragment to pick)
            % EACH LOOP HERE BUILDS A NEW FRAG. GLYPEP
            tempasmbplan = asmbplan(k,:);
            if ~any(tempasmbplan == 0)
                if buildfragsgp
                    newg = g;  % vptm sites in pep. backbone
                end
                newsgpmass = basemz;
                newunitindex = baseunitindex;
                newunitindex_glycan = newunitindex{2};
                newstublen = basestublen;
                iontypeaddon = '';
                for n = 1:length(tempasmbplan)
                    % n MARKS EACH PTM
                    % tempasmbplan(n) MARKS WHICH FRAGMENT TO GRAB &
                    % COMBINE WITH PEP
                    nthtempptmfrag = tempptm_frag{n};
                    nthtempptm_frag_glycanind = tempptm_frag_glycanind{n};
                    if tempptm_typ(n) == 1  % glycan
                        thisptm_frag = nthtempptmfrag(tempasmbplan(n));
                        if buildfragsgp
                            newg(n).struct = thisptm_frag.sgp;
                        end
                        newsgpmass = newsgpmass + thisptm_frag.mz - 18.0105647 - 1.007825032;  % - water - H (because both mz are M+H)
                        newunitindex_glycan = [newunitindex_glycan;...
                            nthtempptm_frag_glycanind{tempasmbplan(n)}];
                            newstublen = [newstublen,thisptm_frag.stublen];
                        if ~strcmpi(thisptm_frag.type,'none')
                            fragtyptxt = thisptm_frag.type;
                            fragtyptxt = strrep(fragtyptxt,'-y-','');
                            fragtyptxt = strrep(fragtyptxt,'none-','');
                            iontypeaddon = [iontypeaddon,fragtyptxt,'-'];
                        end
                    elseif tempptm_typ(n) == 2  % non-glycan
                        if tempptm_frag{n} == 1  % non-glycan mod IS here
                            if buildfragsgp
                                newg(n).struct = tempptm_seq{n};
                            end
                            newsgpmass = newsgpmass + tempptm_mass(n);
                            newunitindex{3} = [newunitindex{3},modunitindex(n)];
                        elseif tempptm_frag{n} == 0  % non-glycan mod is NOT here
                            if buildfragsgp
                                newg(n).struct = '';
                            end
                            iontypeaddon = [iontypeaddon,tempptm_seq{n},'-'];
                        end
                    end
                end
                if ~isempty(iontypeaddon)
                    iontypeaddon = iontypeaddon(1:end-1);
                end
                if strcmpi(basetype,'none')  % intact peptide backbone
                    if ~isempty(iontypeaddon)
                        iontype = ['-y-',iontypeaddon];
                    else
                        iontype = basetype;
                    end
                else
                    if ~isempty(iontypeaddon)
                        iontype = [basetype,'-',iontypeaddon];
                    else
                        iontype = basetype;
                    end
                end
                original{writeind} = baseoriginal;
                if buildfragsgp
                    sgp{writeind} = joinGlyPep(p,newg,m);
                end
                npFrag{writeind} = basenpFrag;
                ngFrag{writeind} = ngnmfragnum(j,1);
                nmFrag{writeind} = ngnmfragnum(j,2);
                mz{writeind} = newsgpmass;
                type{writeind} = iontype;
                charge{writeind} = 1;
                newunitindex{2} = newunitindex_glycan;
                unitindex{writeind} = newunitindex;
                stublen{writeind} = newstublen;
                writeind = writeind + 1;
            end
        end
    end
end
% erase space unused
emptyspaces = false;
if writeind <= estifraglen
    emptyspaces = true;
end
if emptyspaces
    keepind = false(1,length(original));
    keepind(1:writeind-1) = true;
    fragcell = [original(keepind),sgp(keepind),npFrag(keepind),...
        ngFrag(keepind),nmFrag(keepind),mz(keepind),...
        type(keepind),charge(keepind),unitindex(keepind),...
        stublen(keepind)];
else
    fragcell = [original,sgp,npFrag,ngFrag,nmFrag,mz,...
        type,charge,unitindex,stublen];
end
combifrags = cell2struct(fragcell,{'original','sgp','npFrag','ngFrag','nmFrag',...
    'mz','type','charge','unitindex','stublen'},2)';
end