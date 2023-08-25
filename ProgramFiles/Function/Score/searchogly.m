function [result,newprotbatch] = searchogly(thisprotbatch,pepcomp,pepmass,FASTAheadall,...
    pepfrag_HCD,pepfrag_others,fragmode_others,...
    PTM,glycancombi,glycancombimass,scoreoptions,SASSO,colnames,...
    agpoptions,scannum,spectra,charge,precmz,retime,exptmass,precmass,quant,...
    SASSOind)

% This file organizes data for analysis using 'oglypepsearch1spec.m' for
% individual scans

% organize SASSO data for _trigger mode
searchSASSO = SASSO(SASSOind,:);
SASSOcolnames = [colnames,'MS1','ParentRetime'];
HCDcolind = ismember(upper(SASSOcolnames),'HCD');
CIDcolind = ismember(upper(SASSOcolnames),'CID');
EThcDcolind = ismember(upper(SASSOcolnames),'ETHCD');

calcCID = false;
calcEThcD = false;
if any(EThcDcolind)
    calcEThcD = true;
    pepfrags_EThcD = pepfrag_others(:,strcmpi(fragmode_others,'ETHCD'));
end
if any(CIDcolind)
    calcCID = true;
end

% only store a subset of input data corresponding to files in SASSO,
% including MS1 scans
allsearchscannum = searchSASSO(:);
scannumkeepind = ismember(scannum,allsearchscannum);
scannum = scannum(scannumkeepind);
spectra = spectra(scannumkeepind);
charge = charge(scannumkeepind);
precmz = precmz(scannumkeepind);
retime = retime(scannumkeepind);
exptmass = exptmass(scannumkeepind);
precmass = precmass(scannumkeepind);
quant = quant(scannumkeepind);

blocksize = 100000;
fldnms = ResultItem.allresultitems;   % states which fields are output and in what sequence
result = [];
pepseqcounter = 1;
currentind = blocksize + 1;
tempaccuresult = [];
newprotbatch = {};
tempaccuprotbatch = {};

% Now for each SASSO row get HCDspectrum, EThcDspectrum and CIDspectrum
% ALso, collect other data like charge, precms etc. that is needed for
% scoring

for ii = 1:length(searchSASSO(:,1))
%for ii = 1737
    CIDspectrum = {};
    EThcDspectrum = {};
    HCD_scanind = scannum == searchSASSO(ii,HCDcolind);
    HCDspectrum = spectra{HCD_scanind};
    if searchSASSO(ii,HCDcolind) == 4541             % debug step: used to stop for debug at specfic HCD scan
        a=1;
    end
     if calcEThcD                                    % This reads the SASSO file and gets EThcD spectrum
        EThcD_scanind = scannum == searchSASSO(ii,EThcDcolind);
        EThcDspectrum = spectra{EThcD_scanind};
    end
    if calcCID                                      % This reads the SASSO file and gets CID spectrum
        CID_scanind = scannum == searchSASSO(ii,CIDcolind);
        CIDspectrum = spectra{CID_scanind};
    end
    generalcharge = charge(HCD_scanind);
    generalprecmz = precmz(HCD_scanind);
    generalretime = retime(HCD_scanind);
    generalexptmass = exptmass(HCD_scanind);
    generalprecmass = precmass(HCD_scanind);
    generalquant = quant(HCD_scanind);

    % First score the HCD spectrum in workingmode = 1-- results are
    % HCDresults, provided coverage is more than minimum pepcov
    [HCDresult,~] = oglypepsearch1spec(HCDspectrum,searchSASSO(ii,HCDcolind),...
        generalcharge,'HCD',generalprecmz,...
        thisprotbatch,[],pepfrag_HCD,pepcomp,PTM,[],1,scoreoptions); % send the file in working_mode=1 to see if HCD match is sufficient.
    bestHCDglypeprealind = 1:length(thisprotbatch);
    [bestHCDglypepscore,bestHCDglypepind] = sort(HCDresult,'descend');
    bestHCDglypepindmask = false(size(bestHCDglypepind));
    bestHCDglypepindmask(1:min(length(bestHCDglypepind),scoreoptions.keeptopn)) = true;
    bestHCDglypeprealind = bestHCDglypeprealind(bestHCDglypepind);
    bestHCDglypepind = (bestHCDglypepscore >= 0) & bestHCDglypepindmask;
    bestHCDglypeprealind = bestHCDglypeprealind(bestHCDglypepind);
    bestHCDglypepthisbatch = thisprotbatch(bestHCDglypeprealind);
    bestHCDglypepthisbatchmass = pepmass(bestHCDglypeprealind);
    bestHCDglypepthisbatchpepcomp = pepcomp(bestHCDglypeprealind,:);

    % This part finds the combination of peptide and glycans that matches
    % the candidate MS1 mass [precmass(HCD_scanind)]. The result that
    % emerges is the theoretical glypepmass ['tempglypepmass'] using the
    % combination of peptide and glycan described in tempglypepcomp. The
    % first element of each array in the cell array is the peptide number, 
    % the remaining are the modifications (glycan and non-glycan,
    % both fixed and variable) that together can result in the MS1 mass
    % within said tolerance
    [tempglypepcomp,tempglypepmass,glycancombi] = getoglyms1match(...
        bestHCDglypepthisbatch,bestHCDglypepthisbatchmass,...
        precmass(HCD_scanind),PTM,glycancombi,glycancombimass,scoreoptions);

    tempglypepcomp_keepind = true(size(tempglypepcomp));
    for jj = 1:length(tempglypepcomp)
        thistempglypepcomp = tempglypepcomp{jj};
        thistempglypepcomp_PTM = thistempglypepcomp(2:end);
        thistempglypepcomp_PTMisglycan = false(size(thistempglypepcomp_PTM));
        for kk = 1:length(thistempglypepcomp_PTM)
            thistempglypepcomp_PTMisglycan(kk) = PTM.ptmisvar(thistempglypepcomp_PTM(kk));
        end
        uniPTMisglycan = unique(thistempglypepcomp_PTMisglycan);
        if length(uniPTMisglycan) == 1 && ~unique(thistempglypepcomp_PTMisglycan)
            tempglypepcomp_keepind(jj) = false;
        end
    end
    tempglypepcomp = tempglypepcomp(tempglypepcomp_keepind);
    
    if ~isempty(tempglypepcomp)             % if a glypep combination is found localize it
        candipepind = cellfun(@(x) x(1),tempglypepcomp);
        candiPTMind = cellfun(@(x) x(2:end),tempglypepcomp,'Uniformoutput',false);
        candipepseqseq = bestHCDglypepthisbatch(candipepind);
        candipepseqcomp = bestHCDglypepthisbatchpepcomp(candipepind,:);
        if scoreoptions.isHCDtrigger && calcEThcD
            bestHCDthisprotbatchpepfrags_EThcD = pepfrags_EThcD(bestHCDglypeprealind);
            candipepfrag = bestHCDthisprotbatchpepfrags_EThcD(candipepind);
            [thisEThcDresult,thisEThcDresult_comp] = oglypepsearch1spec(EThcDspectrum,searchSASSO(ii,EThcDcolind),generalcharge,...
                'EThcD',generalprecmz,candipepseqseq,tempglypepmass,candipepfrag,candipepseqcomp,PTM,candiPTMind,...
                2,scoreoptions);
            bestpepseq = thisEThcDresult.bestpepseq;
            bestpepmass = thisEThcDresult.bestpepmass;
            bestpepptmcomp = thisEThcDresult.bestpepptmcomp;
            bestpepfragsto = thisEThcDresult.bestpepfragsto;
            bestpepcomp = thisEThcDresult_comp;
        else
            candipepfrag = pepfrag_HCD(candipepind);
            [thisHCDresult,thisHCDresult_comp] = oglypepsearch1spec(HCDspectrum,searchSASSO(ii,HCDcolind),generalcharge,...
                'HCD',generalprecmz,candipepseqseq,tempglypepmass,candipepfrag,PTM,candiPTMind,...
                2,scoreoptions);
            bestpepseq = thisHCDresult.bestpepseq;
            bestpepmass = thisHCDresult.bestpepmass;
            bestpepptmcomp = thisHCDresult.bestpepcomp;
            bestpepcomp = thisHCDresult_comp;
        end

        % Determine the localization type
        fldnms=[fldnms,'localization'];
        if length(bestpepseq) > 0
            nSite=length(strfind(bestpepseq{1},'T'))+length(strfind(bestpepseq{1},'S')); % number of possible O-glycosylation sites
            gPos=[];
            for kk=1:length(bestpepseq)
                [~,g,~]=breakGlyPep(bestpepseq{kk});
                gPos = [gPos, g.pos];
            end
            occupiedSites = length(unique(gPos));
            fracOccupied = occupiedSites/nSite;
            if length(bestpepseq) == 1
                localization = 'localized';
            elseif fracOccupied ==1
                localization = 'Non-localized';
            else
                localization = 'Part-localized';
            end
        end
        % Do scoring for all fragmentation modes
        spectragrp = {HCDspectrum;CIDspectrum;EThcDspectrum};
        spectragrp = spectragrp(ismember({'HCD','CID','ETHCD'},upper(colnames)));
        for jj = 1:length(bestpepseq)
            tempbestpepptmcomp = bestpepptmcomp{jj};
            tempbestpepmass = bestpepmass(jj);
            tempresult = analyzegp({[1,1,tempbestpepptmcomp]},tempbestpepmass,...
                searchSASSO(ii,1:length(colnames)),generalcharge*ones(1,length(colnames)),spectragrp,...
                generalretime*ones(1,length(colnames)),generalexptmass*ones(1,length(colnames)),generalprecmass*ones(1,length(colnames)),...
                generalquant*ones(1,length(colnames)),1:length(colnames),1,{bestpepseq(jj)},FASTAheadall(bestpepcomp(jj,1)),1,PTM,...
                agpoptions);
            %% check in
                for i=1:3
                    tempresult(i).localization = localization;
                end
            %% check out
            for kk = 1:length(tempresult)
                tempProteinID = str2num(tempresult(kk).ProteinID);
                tempProteinID(1) = bestpepcomp(jj,1);
                tempProteinID(2) = pepseqcounter;
                tempresult(kk).ProteinID = num2str(tempProteinID);
                pepseqcounter = pepseqcounter + 1;
                currentind = currentind - 1;
                if currentind == 0
                    result = [result,tempaccuresult];
                    currentind = blocksize;
                    tempaccuresult = [];
                    newprotbatch = [newprotbatch;tempaccuprotbatch];
                    tempaccuprotbatch = {};
                end
                for mm = 1:length(fldnms)
                    tempaccuresult(currentind).(fldnms{mm}) = tempresult(kk).(fldnms{mm});
                end
                tempaccuprotbatch{currentind,1} = bestpepseq{jj};
            end
        end

    end
end
if ~isempty(tempaccuresult)
    tempaccuresult = tempaccuresult(~cellfun(@isempty,{tempaccuresult.Scan}));
    tempaccuprotbatch = tempaccuprotbatch(~cellfun(@isempty,tempaccuprotbatch));
end
result = [result,tempaccuresult];
newprotbatch = [newprotbatch;tempaccuprotbatch];
end