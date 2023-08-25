function [theofragsto,gpfragsto,gpfragstoindex,pepfragsto] = gettheofragsto(goodcandidateind,...
    fragmodeind,goodcomp,prot,ptmseq,ptmfragsto,ptmfragstublen,ptmtype,ptmmass,gpfragsto,...
    gpfragstoindex,pepfragsto,analyzefragmode,fragnum,userpeptiontyp,options)

pepstomode = options.pepstomode;
storeglypepfrag = options.gpfragstomode;
monosacislabile = options.monosacislabile(fragmodeind);
simultaneousfrag = options.simultaneousfrag(fragmodeind);

fragopt.stublen = options.maxstublen;  % maximum stub len on peptide frag
% fragopt.HCDstubfrag_override = true;  % RESERVED: HCD ngfrag be at least 2
fragopt.fragmode = '';
fragopt.monosacislabile = false;
fragopt.simultaneousfrag = false;
fragopt.specialcondition = {};
fragopt.mode = options.sgpseqmode;  % this mode generates(1) SGP for fragment (or not(2))
                    
if storeglypepfrag == 1  
    % store previously generated fragments, reuse if possible VERY HIGH MEMORY USE
    theofragsto = gpfragsto(goodcandidateind,fragmodeind);
    nofrag = cellfun(@isempty,theofragsto);
elseif storeglypepfrag == 0  % rebuild fragments every time
    theofragsto = cell(length(goodcomp),sum(fragmodeind));
end
fragmodeind2 = find(fragmodeind);
for i = 1:length(goodcomp)
    thiscomp_theofragsto = theofragsto(i,:);
    if any(cellfun(@isempty,thiscomp_theofragsto))
        thisgoodcomp = goodcomp{i};
        fragnum_todo = fragnum(fragmodeind,:);
        pepiontyp_todo = userpeptiontyp(fragmodeind,:);
        fragmode_todo = analyzefragmode(fragmodeind);
        for j = 1:sum(fragmodeind)
            if isempty(thiscomp_theofragsto{j})
                %% A. NECESSARY FRAGMENTATION PARAMETER ADJUSTMENTS
                specialcondition = {};  % RESERVED: record what's special for (maybe) special treatment
                thisfragmode = fragmode_todo{j};
                thisfragnum = fragnum_todo(j,:);  % [np ng nm]
                thisfragpepiontyp = pepiontyp_todo{j};
                if strcmpi(thisfragmode,'CID')  % CID special treatment
                    if length(thisgoodcomp) > 2  % with glycan
                        tempPTMseqs = ptmseq(thisgoodcomp(3:end));
                        totalmonosac = sum(cellfun(@(x) length(strfind(x,'{')),tempPTMseqs));
                        if totalmonosac < 6
                            thisfragnum(1) = 1;  % CID, number of monosac too small -> NpFrag = 1
                            specialcondition = [specialcondition;'CID_fragpep,few_monosac'];
                        end
                    else  % without glycan -> NpFrag = 1
                        thisfragnum(1) = 1;
                        specialcondition = [specialcondition;'CID_fragpep,no_monosac'];
                    end
                elseif strcmpi(thisfragmode,'ETCID')
                    if length(thisgoodcomp) == 2  % ETCID without glycan
                        specialcondition = [specialcondition;'ETCID_pepiontyp,no_monosac'];
                    end
                end
                %% B. PREPARE PEPTIDE FRAGMENTS
                pepseq = prot{thisgoodcomp(1)}{thisgoodcomp(2)};
                switch pepstomode
                    case 0  % RE-GENERATE EVERYTIME
                        theofrag_pep = fragpep(pepseq,thisfragnum,thisfragpepiontyp,...
                            thisfragmode);  % frag pep
                    case 1  % EACH FRAGMODE HAS PLACES
                        fragmodeind_todo = find(fragmodeind);
                        theofrag_pep = pepfragsto{fragmodeind_todo(j)}{thisgoodcomp(1)}{thisgoodcomp(2)};
                        if isempty(theofrag_pep)  % build pep frag
                            theofrag_pep = fragpep(pepseq,thisfragnum,thisfragpepiontyp,...
                                thisfragmode);  % frag pep
                            pepfragsto{fragmodeind_todo(j)}{thisgoodcomp(1)}{thisgoodcomp(2)} = theofrag_pep;
                        end
                    case 2  % ONCE, THEN ABRIDGE
                        theofrag_pep = pepfragsto{thisgoodcomp(1)}{thisgoodcomp(2)};
                        if isempty(theofrag_pep)  % build pep frag
                            theofrag_pep = fragpep(pepseq,thisfragnum,'bcyzi',...
                                thisfragmode);  % frag pep
                            pepfragsto{thisgoodcomp(1)}{thisgoodcomp(2)} = theofrag_pep;
                        end
                        theofragpeptyp = cellfun(@(x) x(1),{theofrag_pep.type});  % get ion type
                        if thisfragnum(1) > 0
                            keeplist = arrayfun(@(x) ismember(x,thisfragpepiontyp),theofragpeptyp);
                        else
                            keeplist = false(size(theofragpeptyp));
                        end
                        keeplist(1) = true;  % keep the intact ion, regardless of ion type selection
                        theofrag_pep = theofrag_pep(keeplist);
                end
                if length(thisgoodcomp) > 2  % peptide w/ PTM
                    fragopt.fragmode = thisfragmode;
                    fragopt.monosacislabile = monosacislabile(j);
                    fragopt.simultaneousfrag = simultaneousfrag(j);
                    fragopt.specialcondition = specialcondition;
                    fragopt.addoxoniumion = false;
                    if strcmpi(thisfragmode,'HCD') || strcmpi(thisfragmode,'ETHCD')
                        fragopt.addoxoniumion = true;
                    end
                    temptheofrag = combitheofrag(theofrag_pep,ptmfragsto(fragmodeind2(j),:),...
                        ptmfragstublen(fragmodeind2(j),:),pepseq,ptmseq,...
                        ptmtype,ptmmass,thisgoodcomp,fragnum_todo(j,:),fragopt);
                    [~,ind] = unique(round([temptheofrag.mz],9));
                    temptheofrag = temptheofrag(fliplr(ind));
                    theofragsto{i,j} = temptheofrag;
                else  % peptide w/o PTM
                    [~,ind] = unique([theofrag_pep.mz]);
                    theofrag_pep = theofrag_pep(fliplr(ind));
                    theofragsto{i,j} = theofrag_pep;
                end
            end
        end
    end
end
if storeglypepfrag
    if any(nofrag)
        fragmodeind2 = find(fragmodeind);
        for i = 1:length(fragmodeind2)
            currentgpfragstoindex = gpfragstoindex(:,fragmodeind2(i));
            tempgpfragstoindex = currentgpfragstoindex;
            tempgpfragstoindex(goodcandidateind,:) = true;
            temptheofragsto = theofragsto(:,i);
            gpfragsto(xor(currentgpfragstoindex,tempgpfragstoindex),fragmodeind2(i)) = ...
                temptheofragsto(nofrag(:,i));
            gpfragstoindex(goodcandidateind,fragmodeind2(i)) = true;
        end
    end
end
end