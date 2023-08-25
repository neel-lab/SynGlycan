function [seq,mass,massfix,type,fragsto,stublen] = theoptmfrag(ptmstruct,ngfrag,...
    allfragmode,options)
% THEOPTMFRAG: fragment all PTM appeared in digested protein file
%
% Syntax:
% [seq,mass,massfix,type,fragsto,stublen] =
%     theoptmfrag(ptmstruct,ngfrag,allfragmode,options)
%
% Input:
% (Below "n" equals to number of fragmentation methods.)
% ptmstruct: m x 1 cell array of strings. PTM sequences.
% ngfrag: n x 1 numerical array. Number of cleavages allowed when
%     fragmenting glycans.
% allfragmode: n x 1 cell array of strings. All fragmentation modes that
%     needs to fragment these PTMs for.
% opitons: structure, field(s) are:
%     monosacislabile: n x 1 logical array. If true, labile monosac. "s" and
%     "f" will not be restricted by NgFrag.
%
% Output:
% (Below "n" equals to number of fragmentation methods.)
% seq: Sequence of PTMs.
% mass: Mass of PTMs.
% massfix: Mass adjustments. When trying to get glycopeptide masses, after
%     adding peptide and PTM, add these values to compensate water loss or
%     others.
% type: 1 is glycan, 2 is non-glycan PTM.
% fragsto: n x m cell array of structure. m equals to size of input
%     "ptmstruct". Each cell contains a 1 x k structure of frag. ions.
%     For frag. ion data structure format see MULTISGPFRAG.
%     If PTM is a non-glycan one, the corresponding cells will
%     be empty.
% stublen: n x m cell array. Size is same as "fragsto". The length of
%     glycan fragments. It is necessary for HCD fragmentation.
%
% Note:
% Input does not contain nmFrag (number of non-glycan PTM fragmentation
%     allowed). This is because in current design, fragmentation of non-glycan
%     PTM is done by completely remove this PTM from pepside (a yes/no style).
%     nmFrag is used in later stages of calculation.
% Fragment ions with labile monosac. removed has special marker in field
% Besides "labile monosac.", HCD has special procedure: if NgFrag <= 1 but
%     > 0, glycan Y-ions will be generated using NgFrag = 2.
%
% Example:
% Load file "testdata_theoptmfrag.mat" and run command:
% [seq,mass,massfix,type,fragsto,stublen] = theoptmfrag(ptmstruct,ngfrag,...
%     allfragmode,options)
% to check results
%
% Children function:
% N/A
%
% See also:
% MULTISGPFRAG  GLYMW  PTM  GETDISTANCE

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 04/01/2020

seq = reshape(ptmstruct,1,[]);
mass = zeros(1,length(ptmstruct));
massfix = zeros(1,length(ptmstruct));
type = zeros(1,length(ptmstruct));
fragsto = cell(length(ngfrag),length(ptmstruct));
stublen = cell(size(fragsto));

ngfragfix = ngfrag;
[uningfrag,~,ind] = unique(ngfragfix,'stable');
unifragsto = cell(length(uningfrag),length(ptmstruct));
% if different frag mode have same frag number, then there's no need to do
% it twice, just copy the data
for ii = 1:length(ptmstruct)
    thisptm = ptmstruct{ii};
    if isempty(ptmstruct{ii})
        % INTENTIONALLY LEFT BLANK - THIS HAPPENS WHEN NO PTM IS PRESENT
    elseif strcmp(ptmstruct{ii}(1),'{')  % glycan
        type(ii) = 1;
        for jj = 1:length(uningfrag)
            tempfrags = multiSGPFrag(thisptm,0,uningfrag(jj),0,1,'unimass1');
            % CAUTION: UNIMASS1 IS NOT GOOD FOR MULTIPLE GLYCOSYLATED
            % GLYPEP - LOSE Y-IONS
            [~,tempind] = sort([tempfrags.ngFrag],'descend');
            tempfrags = tempfrags(tempind);
            unifragsto{jj,ii} = tempfrags;
        end
        mass(ii) = glyMW(thisptm);  % glycan mass
        massfix(ii) = -18.0105647;  % water loss when combined
    elseif strcmp(ptmstruct{ii}(1),'<')  % non-glycan, no fragmentation in current version
        type(ii) = 2;
        tempvar = thisptm(2:end-1);
        mass(ii) = ptm(tempvar);  % so far, don't consider mass change during combination
    end
end
% FIND GLYCAN RESIDUE SIZE, THIS IS STUBLEN VALUES FOR SIMULTANEOUS
% FRAGMENTATION
unistublen = cell(size(unifragsto));
for ii = 1:size(unifragsto,1)
    for jj = 1:size(unifragsto,2)
        if ~isempty(unifragsto{ii,jj})
            tempstublen = zeros(size(unifragsto{ii,jj}));
            for k = 1:length(unifragsto{ii,jj})
                tempstublen(k) = max(getdistance(unifragsto{ii,jj}(k).sgp));
            end
            unistublen{ii,jj} = tempstublen;
        end
    end
end

% IN CASE DIFFERENT FRAG MODE HAVE IDENTICAL FRAG NUMBER - COPY THE VALUES
for ii = 1:length(ind)
    fragsto(ii,:) = unifragsto(ind(ii),:);
    stublen(ii,:) = unistublen(ind(ii),:);
end

% HCD special treatment: if ngfrag = 1, remove those B-ions with ngFrag too
% big (NOT IN USE - RESERVED)
% if any(HCDind) && ngfrag(HCDind) == 1
%     for i = 1:length(ptmstruct)
%         tempfrag = fragsto{HCDind,i};
%         tempstublen = stublen{HCDind,i};
%         if ~isempty(tempfrag)
%             alltype_b = cellfun(@(x) strcmpi(x(1:2),'-b'),{tempfrag.type});
%             alltype_none = cellfun(@(x) strcmpi(x,'none'),{tempfrag.type});
%             fragsto{HCDind,i} = [tempfrag(alltype_none),...
%                 tempfrag(([tempfrag.ngFrag] == 1) & alltype_b),...
%                 tempfrag(~(alltype_b | alltype_none))];
%             stublen{HCDind,i} = [tempstublen(alltype_none),...
%                 tempstublen(([tempfrag.ngFrag] == 1) & alltype_b),...
%                 tempstublen(~(alltype_b | alltype_none))];
%         end
%     end
% end

%% GENERATE LABILE FRAGMENTS - REMOVE LABILE MONOSAC 1 BY 1 THEN REDO FRAGMENTATION
% NOT A VERY GOOD DESIGN - BUT WORKS - THERE MAY BE REDUNDANT CALCULATIONS
% HERE CAUSED BY "SHARE NGFRAG"
for ii = 1:length(options.monosacislabile)
    if options.monosacislabile(ii)
        for jj = 1:length(seq)  % each PTM
            thisptm = seq{jj};
            haslabilems = [strfind(thisptm,'f'),strfind(thisptm,'s')];
            if type(jj) == 1 && ~isempty(haslabilems)
                indtemp = strfind(thisptm,'{');
                levelindex = zeros(2,length(thisptm));
                indtemp2 = strfind(thisptm,'}');
                levelindex(1,indtemp) = 1;
                levelindex(1,indtemp2) = -1;
                levelindex(2,1) = 1;
                for k = 2:size(levelindex,2)
                    levelindex(2,k) = levelindex(2,k-1) + levelindex(1,k);
                end
                haslabilems = sort(haslabilems,'descend');
                msposinseq = find(levelindex(1,:) == 0);
                for k = length(haslabilems):-1:1  % remove 1 ~ max number of fragile MS
                    delmsplan = combnk(haslabilems,k);
                    for m = 1:size(delmsplan,1)  % each removal plan
                        thisdelmsplan = delmsplan(m,:);
                        tempmsposinseq = msposinseq;
                        tempmsposinseq(ismember(tempmsposinseq,thisdelmsplan)) = 0;
                        removedmonosac = find(tempmsposinseq == 0);
                        fixedmsind = setdiff(1:length(msposinseq),removedmonosac);
                        tempptm = thisptm;
                        newlevelindex = levelindex(2,:);
                        labilemslevel = newlevelindex(thisdelmsplan);
                        labilemssgp = tempptm(thisdelmsplan);
                        for n = 1:length(labilemslevel)  % each fMS to be removed
                            startind = labilemslevel(n);
                            tgtind = startind - 1;
                            tgtpos = find(newlevelindex(thisdelmsplan(n):end) == tgtind,1) + thisdelmsplan(n) - 1;
                            tempptm([tgtpos,thisdelmsplan(n)-1,thisdelmsplan(n)]) ='';
                            newlevelindex([tgtpos,thisdelmsplan(n)-1,thisdelmsplan(n)]) =[];
                        end
                        labiledfrags = multiSGPFrag(tempptm,0,ngfragfix(ii),0,1,'allion');
                        [~,ind] = sort([labiledfrags.ngFrag],'descend');
                        labiledfrags = labiledfrags(ind);
                        labiledmonosac = '';
                        for n = 1:length(labilemssgp)
                            labiledmonosac = [labiledmonosac,'-{',labilemssgp(n),'}'];
                        end
                        labiledstublen = zeros(size(labiledfrags));
                        for n = 1:length(labiledfrags)
                            labiledfrags(n).type = [labiledfrags(n).type,'-labile',labiledmonosac];
                            labiledfrags(n).unitindex{2}(:,1) = fixedmsind(labiledfrags(n).unitindex{2}(:,1));
                            labiledstublen(n) = max(getdistance(labiledfrags(n).sgp));
                        end
                        fragsto{ii,jj} = [fragsto{ii,jj},labiledfrags];
                        stublen{ii,jj} = [stublen{ii,jj},labiledstublen];
                    end
                end
                theofrags = fragsto{ii,jj};
                stublens = stublen{ii,jj};
                
                alliontyp = {theofrags.type};
                refunitind_gly = theofrags(1).unitindex{2}(end,1);
                isintact = cellfun(@(x) strcmpi(x(1:4),'none'),alliontyp);
                theofrags_intact = theofrags(isintact);
                stublens_intact = stublens(isintact);
                unitinds_intact = zeros(length(theofrags_intact),refunitind_gly);
                for k = 1:length(theofrags_intact)
                    unitinds_intact(k,theofrags_intact(k).unitindex{2}(:,1)) = 1;
                end
                [~,ind_intact] = unique(unitinds_intact,'rows','stable');
                theofrags_intact = theofrags_intact(ind_intact);
                [theofrags_intact.ngFrag] = deal(1);
                theofrags_intact(1).ngFrag = 0;
                stublens_intact = stublens_intact(ind_intact);
                
                isbion = cellfun(@(x) strcmpi(x(1:2),'-b'),alliontyp);
                theofrags_b = theofrags(isbion);
                stublens_b = stublens(isbion);
                unitinds_b = zeros(length(theofrags_b),refunitind_gly);
                for k = 1:length(theofrags_b)
                    unitinds_b(k,theofrags_b(k).unitindex{2}(:,1)) = 1;
                end
                [~,ind_b] = unique(unitinds_b,'rows','stable');
                theofrags_b = theofrags_b(ind_b);
                stublens_b = stublens_b(ind_b);
                
                isyion = ~(isbion | isintact);
                theofrags_y = theofrags(isyion);
                [~,indngf_y] = sort([theofrags_y.ngFrag],'descend');
                theofrags_y = theofrags_y(indngf_y);
                stublens_y = stublens(isyion);
                stublens_y = stublens_y(indngf_y);
                unitinds_y = zeros(length(theofrags_y),refunitind_gly);
                for k = 1:length(theofrags_y)
                    unitinds_y(k,theofrags_y(k).unitindex{2}(:,1)) = 1;
                end
                [~,ind_y] = unique(unitinds_y,'rows','stable');
                theofrags_y = theofrags_y(ind_y);
                stublens_y = stublens_y(ind_y);
                
                newtheofrags = [theofrags_intact,theofrags_b,theofrags_y];
                newstublens = [stublens_intact,stublens_b,stublens_y];
                fragsto{ii,jj} = newtheofrags;
                stublen{ii,jj} = newstublens;
            end
        end
    end
    
end
end