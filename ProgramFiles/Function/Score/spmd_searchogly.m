function [result,newprotbatch] = spmd_searchogly(thisprotbatch,thisprotbatch_pepfrag,...
    FASTAbatch,thisprotbatchmass,PTM,glycancombi,glycancombimass,scoreoptions,SASSO,colnames,...
    agpoptions,scannum,spectra,charge,precmz,retime,exptmass,precmass,quant,dist_SASSOind)

searchSASSO = SASSO(dist_SASSOind,:);
SASSOcolnames = [colnames,'MS1','ParentRetime'];
HCDcolind = ismember(upper(SASSOcolnames),'HCD');
CIDcolind = ismember(upper(SASSOcolnames),'CID');
EThcDcolind = ismember(upper(SASSOcolnames),'ETHCD');
newprotbatch = {};
result = [];
pepseqcounter = 1;
for ii = 1:length(searchSASSO(:,1))
    HCD_scanind = scannum == searchSASSO(ii,HCDcolind);
    CID_scanind = scannum == searchSASSO(ii,CIDcolind);
    EThcD_scanind = scannum == searchSASSO(ii,EThcDcolind);
    HCDspectrum = {};
    CIDspectrum = {};
    EThcDspectrum = {};
    if ~isempty(HCD_scanind)
        HCDspectrum = spectra{HCD_scanind};
    end
    if ~isempty(EThcD_scanind)
        EThcDspectrum = spectra{EThcD_scanind};
    end
    if ~isempty(CID_scanind)
        CIDspectrum = spectra{CID_scanind};
    end
    generalcharge = charge(HCD_scanind);
    generalprecmz = precmz(HCD_scanind);
    generalretime = retime(HCD_scanind);
    generalexptmass = exptmass(HCD_scanind);
    generalprecmass = precmass(HCD_scanind);
    generalquant = quant(HCD_scanind);
    thisprotbatchpepfrags_HCD = cell(size(thisprotbatch));
    thisprotbatchpepfrags_EThcD = cell(size(thisprotbatch));
    for jj = 1:length(thisprotbatch)
        temppepfrag = thisprotbatch_pepfrag{jj};
        thisprotbatchpepfrags_HCD{jj} = temppepfrag(~cellfun(@(x) any(ismember(x,'cz')),...
            {temppepfrag.type}));
        thisprotbatchpepfrags_EThcD{jj} = temppepfrag(~cellfun(@(x) any(ismember(x,'byz')),...
            {temppepfrag.type}));
    end
    [HCDresult,~] = oglypepsearch1spec(HCDspectrum,searchSASSO(ii,HCDcolind),...
        generalcharge,'HCD',generalprecmz,...
        thisprotbatch,[],thisprotbatchpepfrags_HCD,PTM,[],[],1,scoreoptions);
    bestHCDglypeprealind = 1:length(thisprotbatch);
    [bestHCDglypepscore,bestHCDglypepind] = sort(HCDresult,'descend');
    bestHCDglypepindmask = false(size(bestHCDglypepind));
    bestHCDglypepindmask(1:min(length(bestHCDglypepind),scoreoptions.keeptopn)) = true;
    bestHCDglypeprealind = bestHCDglypeprealind(bestHCDglypepind);
    bestHCDglypepind = (bestHCDglypepscore >= 0) & bestHCDglypepindmask;
    bestHCDglypeprealind = bestHCDglypeprealind(bestHCDglypepind);
    
    bestHCDglypepthisbatch = thisprotbatch(bestHCDglypeprealind);
    bestHCDFASTAthisbatch = FASTAbatch;
    bestHCDglypepthisbatchmass = thisprotbatchmass(bestHCDglypeprealind);
    bestHCDthisprotbatchpepfrags_EThcD = thisprotbatchpepfrags_EThcD(bestHCDglypeprealind);
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
    
    if ~isempty(tempglypepcomp)
        candipepind = cellfun(@(x) x(1),tempglypepcomp);
        candiPTMind = cellfun(@(x) x(2:end),tempglypepcomp,'Uniformoutput',false);
        candipepseqseq = bestHCDglypepthisbatch(candipepind);
        candipepfrag = bestHCDthisprotbatchpepfrags_EThcD(candipepind);
        bestglypepthisbatch = bestHCDglypepthisbatch;
        if scoreoptions.isHCDtrigger
            if ~isempty(EThcD_scanind)
                [thisEThcDresult,~] = oglypepsearch1spec(EThcDspectrum,searchSASSO(ii,EThcDcolind),generalcharge,...
                    'EThcD',generalprecmz,candipepseqseq,tempglypepmass,candipepfrag,PTM,candiPTMind,...
                    [],2,scoreoptions);
                bestpepseq = thisEThcDresult.bestpepseq;
                bestpepmass = thisEThcDresult.bestpepmass;
                bestpepcomp = thisEThcDresult.bestpepcomp;
            end
        else
            [thisHCDresult,~] = oglypepsearch1spec(HCDspectrum,searchSASSO(ii,HCDcolind),generalcharge,...
                'HCD',generalprecmz,candipepseqseq,tempglypepmass,candipepfrag,PTM,candiPTMind,...
                [],2,scoreoptions);
            bestpepseq = thisHCDresult.bestpepseq;
            bestpepmass = thisHCDresult.bestpepmass;
            bestpepcomp = thisHCDresult.bestpepcomp;
        end
        spectragrp = {HCDspectrum;CIDspectrum;EThcDspectrum};
        spectragrp = spectragrp(ismember({'HCD','CID','ETHCD'},upper(colnames)));
        for jj = 1:length(bestpepseq)
            tempbestpepcomp = bestpepcomp{jj};
            tempbestpepmass = bestpepmass(jj);
            bestglypepthisbatch{tempbestpepcomp(1)} = bestpepseq{jj};
            tempresult = analyzegp({[1,tempbestpepcomp]},tempbestpepmass,...
                searchSASSO(ii,1:length(colnames)),generalcharge*ones(1,length(colnames)),spectragrp,...
                generalretime*ones(1,length(colnames)),generalexptmass*ones(1,length(colnames)),generalprecmass*ones(1,length(colnames)),...
                generalquant*ones(1,length(colnames)),1:length(colnames),1,{bestglypepthisbatch},{bestHCDFASTAthisbatch},1,PTM,...
                agpoptions);
            newprotbatch = [newprotbatch;bestglypepthisbatch{tempbestpepcomp(1)}];
            for mm = 1:length(tempresult)
                tempresult(mm).ProteinID(2) = pepseqcounter;
            end
            pepseqcounter = pepseqcounter + 1;
            result = [result,tempresult];
        end
    end
end
end