function [result,result_pepcomp] = oglypepsearch1spec(spectrum,scannum,charge,fragmode,precmz,...
    pepseq,pepmass,pepfrag,pepcomp,PTM,PTMcomp,workingmode,scoreoptions)
% OGLYPEPSEARCH1SPEC: search for O linked glycopeptide 1 spectrum at a
%     time, multiple peptide input is allowed
% PTMcombi: n x 1 cell array of numeric array, possible glycan combinations
% pepseq: n x 1 cell array
% pepfrag: n x 1 cell array
%
% The outputs are 'thisEThcDresult' and 'thisEThcDresult_comp'

% This program works in two modes:
% workingmode 1 for HCD backbone matching
% working mode 2 for ETHCD spectrum analysis
% In this mode the two key input parameters are: i) pepseq which is the
% candidate peptide and PTMcomp which is the list of correspinding PTMs
% that sum up to result in precursormz

HCDpepcovthreshold = GlycoPATConstant.Pepcovcutoff_Oglysearch;
EThcDpepcovthreshold = GlycoPATConstant.Pepcovcutoff_Oglysearch;

emp = scoreoptions.analyzefragmode;
if isempty(emp) == 1
    disp(scannum)
end

allfragmodes = scoreoptions.analyzefragmode;
allms2tol = scoreoptions.ms2tol;
allms2tolunit = scoreoptions.ms2tolunit;
methodind = ismember(upper(allfragmodes),upper(fragmode));
ms2tol = allms2tol(methodind);
ms2tolunit = allms2tolunit{methodind};
PTMuniind = PTM.ptmuniind;
PTMisvar = PTM.ptmisvar;
PTMisgly = PTM.ptmisglycan;
PTMrealmass = PTM.ptmmass + PTM.ptmmassfix;
if scannum==3269
    aaa=1
end
switch workingmode
    case 1  % PEP FRAG ONLY
        result = zeros(size(pepfrag));
        result_pepcomp = zeros(size(pepfrag,1),2);
        for ii = 1:length(pepfrag)
            thispepfrag = pepfrag{ii};
            thispepfrag = thispepfrag(~strcmpi({thispepfrag.type},'none'));
            [spectrummatchind,fragmatchind] = quickspecmatch(spectrum,...
                [thispepfrag.mz],charge,ms2tol,ms2tolunit);
            pepcov = fragmatchind(2:2:end) | fragmatchind(1:2:end);
            if sum(pepcov) / length(pepcov) >= HCDpepcovthreshold
                result(ii) = sum(spectrum(spectrummatchind,2));
                result_pepcomp(ii,:) = pepcomp(ii,:);
            end
            % result is total intensity of peaks matched
        end
    case 2  % CONSIDER PTM MASS - STEP BY STEP
        precmz1 = (precmz - 1.007825032) * charge + 1.007825032;
        pepscoresto = [];
        pepseqcovsto = [];
        pepseqindsto = [];
        pepmasssto = [];
        pepcompindsto = [];
        ptmlayoutsto = {};
        pepfragsto = {};  % where are the bonds broken
        for ii = 1:length(pepseq)  % EACH CANDIDATE PEP
            thispepseq = pepseq{ii};
            [p,~,~] = breakGlyPep(thispepseq);
            peplen = length(p.pep);
            thisPTMcomp = PTMcomp{ii};
            if ~isempty(thisPTMcomp)
                thispepfrag = pepfrag{ii};
                thispepfrag = thispepfrag(~strcmpi({thispepfrag.type},'none'));
                [pepseq_PTMposstart,pepseq_PTMposend,pepseq_PTMseq] = ...
                    regexp(thispepseq,'\{\d+\}','start','end','match');  % find and locate PTMs
                ptmpos = pepseq_PTMposstart - 1;  % aaind where ptm resides
                ptmuniind_onpep = zeros(size(pepseq_PTMseq));
                for jj = 1:length(pepseq_PTMseq)
                    ptmuniind_onpep(jj) = str2double(pepseq_PTMseq{jj}(2:end-1));  % 1, 2, 3,... represents PTM marker
                end  % not suitable for determining gly/non-gly, use it with PTM.ptmuniind
                if length(ptmpos) > 1
                    ptmseqlen = pepseq_PTMposend - pepseq_PTMposstart + 1;
                    ptmseqlenfix = ptmseqlen(1:end - 1) * triu(ones(length(ptmseqlen) - 1));
                    ptmpos(2:end) = ptmpos(2:end) - ptmseqlenfix;
                end
                ptmpos_isfix = [];
                ptmpos_isvarng = [];
                ptmpos_isvargly = [];
                for jj = 1:length(pepseq_PTMseq)
                    tempptmuniind = ptmuniind_onpep(jj);
                    if any(PTMisvar(PTMuniind == tempptmuniind))
                        if any(PTMisgly(PTMuniind == tempptmuniind))
                            ptmpos_isvargly = [ptmpos_isvargly,ptmpos(jj)];
                        else
                            ptmpos_isvarng = [ptmpos_isvarng,ptmpos(jj)];
                        end
                    else
                        ptmpos_isfix = [ptmpos_isfix,ptmpos(jj)];
                    end
                end
                ptmcomp_isfix = [];
                ptmcomp_isvarng = [];
                ptmcomp_isvargly = [];
                for jj = 1:length(thisPTMcomp)
                    if PTMisvar(thisPTMcomp(jj))
                        if PTMisgly(thisPTMcomp(jj))
                            ptmcomp_isvargly = [ptmcomp_isvargly,thisPTMcomp(jj)];
                        else
                            ptmcomp_isvarng = [ptmcomp_isvarng,thisPTMcomp(jj)];
                        end
                    else
                        ptmcomp_isfix = [ptmcomp_isfix,thisPTMcomp(jj)];
                    end
                end
                ptmlayout = [];
                subptmpos_isvargly = nchoosek(ptmpos_isvargly,length(ptmcomp_isvargly));
                ptmelem_isvargly = unique(perms(ptmcomp_isvargly),'rows');
                if isempty(ptmpos_isvarng)
                    for jj = 1:size(subptmpos_isvargly,1)
                        for kk = 1:size(ptmelem_isvargly,1)
                            tempptmlayout_isvargly = zeros(1,peplen);
                            tempptmlayout_isvargly(subptmpos_isvargly(jj,:)) = ptmelem_isvargly(kk,:);
                            ptmlayout = [ptmlayout;tempptmlayout_isvargly];
                        end
                    end
                else
                    ptmlayout_isvarng = nchoosek(ptmpos_isvarng,length(ptmcomp_isvarng));
                    subptmpos_isvarng = perms(ptmcomp_isvarng);
                    for jj = 1:size(ptmlayout_isvarng,1)
                        for kk = 1:size(subptmpos_isvarng,1)
                            tempptmlayout_isvar = zeros(1,peplen);
                            tempptmlayout_isvar(ptmlayout_isvarng(jj,:)) = subptmpos_isvarng(kk,:);
                            for ll = 1:size(subptmpos_isvargly,1)
                                for mm = 1:size(ptmelem_isvargly,1)
                                    tempptmlayout_isvar2 = tempptmlayout_isvar;
                                    tempptmlayout_isvar2(subptmpos_isvargly(ll,:)) = ptmelem_isvargly(mm,:);
                                    ptmlayout = [ptmlayout;tempptmlayout_isvar2];
                                end
                            end
                        end
                    end
                end
                for jj = 1:size(ptmlayout,1)
                    ptmlayout(jj,ptmpos_isfix) = ptmcomp_isfix;
                end
                temppepseqcov = zeros(size(ptmlayout,1),1);
                temppepintscore = zeros(size(ptmlayout,1),1);
                searchpepfrag = thispepfrag(cellfun(@(x) any(ismember(x,'c')),...
                    {thispepfrag.type}));
                tempmz_forward = [searchpepfrag.mz];
                for jj = 1:size(ptmlayout,1)
                    tempmz_forward_ptm = tempmz_forward;
                    for kk = 1:length(searchpepfrag)
                        tempaaind = searchpepfrag(kk).unitindex{1};
                        tempptmonpepfrag = ptmlayout(jj,tempaaind);
                        tempptmonpepfrag = tempptmonpepfrag(tempptmonpepfrag > 0);
                        tempptmmassaddon = sum(PTMrealmass(tempptmonpepfrag));
                        tempmz_forward_ptm(kk) = tempmz_forward_ptm(kk) + tempptmmassaddon;
                    end
                    [specmatch_forward,fragmatch_forward] = quickspecmatch(spectrum,tempmz_forward_ptm,charge,...
                        ms2tol,ms2tolunit);         % This is to compare the c-ions
                    [specmatch_reverse,fragmatch_reverse] = quickspecmatch(spectrum,precmz1 - tempmz_forward_ptm + 2*1.007825032,...
                        charge,ms2tol,ms2tolunit);  % This is to compare the z-ions
                    temppepintscore(jj) = sum(spectrum(specmatch_forward,2)) + ...
                        sum(spectrum(specmatch_reverse,2));
                    %                     temppepseqcov(jj) =
                    %                     sum(fragmatch_forward |
                    %                     fragmatch_reverse)/length(fragmatch_forward);
                    %                     % THE ABOVE METHOD COMBINES c- AND
                    %                     z-ION 
                    temppepseqcov(jj) = (sum(fragmatch_forward)/length(fragmatch_forward) + ...
                        sum(fragmatch_reverse)/length(fragmatch_reverse))/2;
                    % THE ABOVE METHOD TREAT c- AND z-IONS SEPARATELY
                    ptmlayoutsto = [ptmlayoutsto;{ptmlayout(jj,:)}];
                    pepfragsto = [pepfragsto; int8(or(fragmatch_forward,fragmatch_reverse))];  % list of found bonds
                end
                pepscoresto = [pepscoresto;temppepintscore];
                pepseqcovsto = [pepseqcovsto;temppepseqcov];
                pepseqindsto = [pepseqindsto;repmat(ii,size(ptmlayout,1),1)];
                pepmasssto = [pepmasssto;repmat(pepmass(ii),size(ptmlayout,1),1)];
            else  % if thisPTMcomp is empty
                pepscoresto = -1;
                pepseqcovsto = -1;
                pepseqindsto = [];
                pepmasssto = -1;
            end
        end   % for each candidate glypep

        % keep only data that are above EThcDpepcovthreshold
        goodpepcovind = pepseqcovsto >= EThcDpepcovthreshold;
        goodpepcovind = and((pepseqcovsto >= EThcDpepcovthreshold),pepseqcovsto == max(pepseqcovsto));   % only keep the hits with maximum coverage. This will reduce the size of the 1_1 file and prepare it for bestisomer
        goodpepscoresto = pepscoresto(goodpepcovind);
        goodpepseqcovsto = pepseqcovsto(goodpepcovind);
        goodpepseqindsto = pepseqindsto(goodpepcovind);
        goodptmlayoutsto = ptmlayoutsto(goodpepcovind);
        goodpepmasssto = pepmasssto(goodpepcovind);
        goodpepfragsto = pepfragsto(goodpepcovind);

        [besttemppepindscore,bestpepscoreind] = sort(goodpepscoresto,'descend');
        bestptmcompind = bestpepscoreind(1:min(length(bestpepscoreind),scoreoptions.keeptopn));   % keeptopn is set in GlycoPATConstant and is currently very large
        bestpepseqcovsto = goodpepseqcovsto(bestptmcompind);
        bestpepseqindsto = goodpepseqindsto(bestptmcompind);
        bestptmlayout = goodptmlayoutsto(bestptmcompind);
        bestpepmasssto = goodpepmasssto(bestptmcompind);
        bestpepfragsto = goodpepfragsto(bestptmcompind);

        bestpepseq = cell(size(bestpepseqindsto));
        bestpepptmcomp = cell(size(bestpepseqindsto));
        for ii = 1:length(bestpepseqindsto)
            thispepseq = pepseq{bestpepseqindsto(ii)};
            bestptmcomp = bestptmlayout{ii};
            [p,g,m] = breakGlyPep(thispepseq);
            gpos = [g.pos];
            keepg = false(size(g));
            for jj = 1:length(g)
                if bestptmcomp(gpos(jj)) ~= 0
                    keepg(jj) = true;
                end
            end
            bestpepseq{ii} = joinGlyPep(p,g(keepg),m);
            bestpepptmcomp{ii} = bestptmcomp(bestptmcomp > 0);
        end
        result.bestpepseq = bestpepseq;
        result.bestpepptmcomp = bestpepptmcomp;
        result.bestpepseqcov = bestpepseqcovsto;
        result.bestpepscore = besttemppepindscore;
        result.bestpepmass = bestpepmasssto;
        result.bestpepfragsto = bestpepfragsto;     % returns information on what broken bonds are seen
        result_pepcomp = pepcomp(bestpepseqindsto,:);
end
end

%% SIMPLE VERSION: HANDLES 1 FRAGMZ INPUT ONLY
function fragmatch = quickspecmatch_simp(spectrum,fragmz,charge,ms2tol,ms2tolunit)
allfragmz = (fragmz - 1.007825032)/(1:charge) + 1.007825032;
fragmatch =  false(size(spectrum,1),1);
switch upper(ms2tolunit)
    case 'DA'
        for ii = 1:length(allfragmz)  % each charge state
            peakhit = any(abs(spectrum(:,1) - allfragmz(ii)) <= ms2tol);
            fragmatch = fragmatch | peakhit;
        end
    case 'PPM'
        for ii =  1:length(allfragmz)
            peakhit = any(abs(spectrum(:,1) - allfragmz(ii))/allfragmz(ii) * 1e6 <= ms2tol);
            fragmatch = fragmatch | peakhit;
        end
end
end

%% MORE CMPLX VERSION: HANDLES MULTIPLE FRAGMZ INPUT
function [specmatch,fragmatch] = quickspecmatch(spectrum,fragmz,charge,ms2tol,ms2tolunit)
allfragmz = (fragmz(:) - 1.007825032)./(1:charge) + 1.007825032;
% (ABOVE) ALLFRAGMZ is a m x n numerical array, for element [a,b] in this array, a
%     means ath fragment, b means charge = b;
specmatch =  false(size(spectrum,1),1);
fragmatch = false(length(fragmz),1);
switch upper(ms2tolunit)
    case 'DA'
        for ii =  1:size(allfragmz,1)  % each frag
            for jj = 1:size(allfragmz,2)  % each charge state
                peakhit = abs(spectrum(:,1) - allfragmz(ii,jj)) <= ms2tol;
                if any(peakhit)
                    specmatch = specmatch | peakhit;
                    fragmatch(ii) = true;
                end
            end
        end
    case 'PPM'
        for ii =  1:size(allfragmz,1)
            for jj = 1:size(allfragmz,2)
                peakhit = abs(spectrum(:,1) - allfragmz(ii,jj))/allfragmz(ii,jj) * 1e6 <= ms2tol;
                if any(peakhit)
                    specmatch = specmatch | peakhit;
                    fragmatch(ii) = true;
                end
            end
        end
end
end