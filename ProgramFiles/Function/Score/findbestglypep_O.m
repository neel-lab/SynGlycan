function findbestglypep_O(input,allfragments,agpoptions,scannum,spectra,allpepdata)

%% LOAD RESULT DATA
numprotperbatch = input.doparacomp * input.numprotperworker;
totalprotnum = length(allpepdata.glypep);
totalfilenum = ceil(totalprotnum/numprotperbatch);
[~,f,e] = fileparts(input.outputfname);
p = input.outputdir;
tgtfolderdir = dir(p);
existingdatafilenames = {tgtfolderdir.name};
existingdatafilenames = existingdatafilenames(~[tgtfolderdir.isdir]);
for ii = 1:length(existingdatafilenames)
    existingdatafilenames{ii} = fullfile(input.outputdir,existingdatafilenames{ii});
end
accumulatedresult = [];
for ii = 1:totalfilenum
    protindse = [(ii-1)*numprotperbatch+1,...
        min(ii*numprotperbatch,totalprotnum)];
    tempresultfilename = fullfile(p,[f,'_',num2str(protindse(1)),'_',num2str(protindse(2)),e]);
    if ismember(tempresultfilename,existingdatafilenames)
        if ii == 1
            load(tempresultfilename,'result','scoreintdata');
        else
            load(tempresultfilename,'result');
        end
        accumulatedresult = [accumulatedresult,result];
    else
        tempresultfilename_regexppattern = [f,'_',num2str(protindse(1)),'_',...
            num2str(protindse(2)),'_Part_[0-9]+of[0-9]+',e];
        fileexists = ~cellfun(@isempty,regexp(existingdatafilenames,tempresultfilename_regexppattern));
        temppresearchfilenames = existingdatafilenames(fileexists);
        for jj = 1:length(temppresearchfilenames)
            tempresultfilename = temppresearchfilenames{jj};
            if (ii == 1)  && (jj == 1)
                load(tempresultfilename,'result','scoreintdata');
            else
                load(tempresultfilename,'result');
            end
            accumulatedresult = [accumulatedresult,result];
        end
    end
end
result = accumulatedresult;
clear accumulate*;

%% RUN DECOY GLYPEP AGAINST RESULT
scoreoptions = scoreintdata.scoreoptions;
denoisingoptHCD.cutoffmed = scoreoptions.cutoffmed(...
    ismember(upper(scoreoptions.analyzefragmode),'HCD'));
denoisingoptHCD.fracmax = scoreoptions.fracmax(...
    ismember(upper(scoreoptions.analyzefragmode),'HCD'));
denoisingoptHCD.fragmode = 'HCD';
denoisingoptCID.cutoffmed = scoreoptions.cutoffmed(...
    ismember(upper(scoreoptions.analyzefragmode),'CID'));
denoisingoptCID.fracmax = scoreoptions.fracmax(...
    ismember(upper(scoreoptions.analyzefragmode),'CID'));
denoisingoptCID.fragmode = 'CID';
denoisingoptETD.cutoffmed = scoreoptions.cutoffmed(...
    ismember(upper(scoreoptions.analyzefragmode),'ETD'));
denoisingoptETD.fracmax = scoreoptions.fracmax(...
    ismember(upper(scoreoptions.analyzefragmode),'ETD'));
denoisingoptETD.fragmode = 'ETD';
denoisingoptEThcD.cutoffmed = scoreoptions.cutoffmed(...
    ismember(upper(scoreoptions.analyzefragmode),'ETHCD'));
denoisingoptEThcD.fracmax = scoreoptions.fracmax(...
    ismember(upper(scoreoptions.analyzefragmode),'ETHCD'));
denoisingoptEThcD.fragmode = 'EThcD';
denoisingoptETciD.cutoffmed = scoreoptions.cutoffmed(...
    ismember(upper(scoreoptions.analyzefragmode),'ETCID'));
denoisingoptETciD.fracmax = scoreoptions.fracmax(...
    ismember(upper(scoreoptions.analyzefragmode),'ETCID'));
denoisingoptETciD.fragmode = 'ETciD';
if scoreintdata.scoreoptions.doparallel
    if scoreintdata.scoreoptions.isHCDtrigger
        protids = {result.ProteinID};
        scans = [result.Scan];
        displaydataind = grouptriggeredresult(scoreintdata.SASSO,protids,scans,'regular');
    else
        displaydataind = (1:length(result))';
    end
    resultind = 1:size(displaydataind,1);
    dist_resultind_dist = distributed(resultind);
    spmd
        dist_result_rev = recalcrevresult(result,displaydataind,allfragments,getLocalPart(dist_resultind_dist),...
            spectra,scannum,...
            denoisingoptCID,denoisingoptHCD,denoisingoptETD,denoisingoptETciD,denoisingoptEThcD,...
            scoreintdata,agpoptions);
    end
    result_rev = [];
    for ii = 1:length(dist_result_rev)
        result_rev = [result_rev,dist_result_rev{ii}];
    end
else
    if scoreintdata.scoreoptions.isHCDtrigger
        protids = {result.ProteinID};
        scans = [result.Scan];
        displaydataind = grouptriggeredresult(scoreintdata.SASSO,protids,scans,'regular');
    else
        displaydataind = (1:length(result))';
    end
    resultind = 1:size(displaydataind,1);
    result_rev = recalcrevresult(result,displaydataind,allfragments,resultind,spectra,scannum,...
        denoisingoptCID,denoisingoptHCD,denoisingoptETD,denoisingoptETciD,denoisingoptEThcD,...
        scoreintdata,agpoptions);
end
protids = {result.ProteinID};
scans = [result.Scan];
if scoreintdata.sliminput.dohcdtrigger
    displaydataind = grouptriggeredresult(scoreintdata.SASSO,protids,scans,'regular');
else
    displaydataind =  (1:length(result))';
end
[CombiES,~] = calculateCombiES(result,displaydataind,...
    scoreintdata.colnames);
protids = {result_rev.ProteinID};
scans = [result_rev.Scan];
if scoreintdata.sliminput.dohcdtrigger
    displaydataind = grouptriggeredresult(scoreintdata.SASSO,protids,scans,'regular');
else
    displaydataind =  (1:length(result_rev))';
end
[CombiES_rev,~] = calculateCombiES(result_rev,displaydataind,...
    scoreintdata.colnames);

enscorecounter = zeros(101,1);
fdrenscorecounter = zeros(101,1);
for ii = 1:101
    enscorecounter(ii) = sum(CombiES >= (ii-1)*.01);
    fdrenscorecounter(ii) = sum(CombiES_rev >= (ii-1)*.01);
end
FDRratio = fdrenscorecounter ./ enscorecounter;
FDRratio(arrayfun(@isnan,FDRratio)) = 0;
badFDR = FDRratio > GlycoPATConstant.OGlyPepFDR;
% result from FIND still fails FDR cut off but + 1 will get you a good one
CombiEScutoff = ones(1,size(FDRratio,2));
for ii = 1:size(FDRratio,2)
    CombiEScutoff(ii) = (find(badFDR(:,ii),1,'last') + 1) / 100;
end
CombiEScutoff(CombiEScutoff > 1) = 1;
displaydataind_keepind = CombiES >= CombiEScutoff;
displaydataind_filtered = displaydataind(displaydataind_keepind,:);
result_filtered = result(displaydataind_filtered(:));
protids = {result_filtered.ProteinID};
scans = [result_filtered.Scan];
displaydataind_fdr = grouptriggeredresult(scoreintdata.SASSO,protids,scans,'COMBINEISOMER');
[CombiES_fdr,~] = calculateCombiES(result_filtered,displaydataind_fdr,...
    scoreintdata.colnames);
displaydataind_postcluster = displaydataind_fdr;
for ii = 1:size(displaydataind_fdr,1)
    tempCombiES_fdr = CombiES_fdr{ii};
    if length(tempCombiES_fdr) > 1 && ...
            max(tempCombiES_fdr) - min(tempCombiES_fdr) > 0.1
        [idx,~] = kmeans(CombiES_fdr{ii}(:),2);
        [uniidx,~,uniidxind] = unique(idx);
        uniidx_meanscore = zeros(size(uniidx));
        for jj = 1:length(uniidx)
            uniidx_meanscore(jj) = mean(CombiES_fdr{ii}(uniidxind == jj));
        end
        [~,highscoreind] = max(uniidx_meanscore);
        tgtidx = uniidx(highscoreind);
        for jj = 1:size(displaydataind_fdr,2)
            displaydataind_postcluster{ii,jj} = displaydataind_postcluster{ii,jj}(idx == tgtidx);
        end
    end
end
displaydataind_final_cell = displaydataind_postcluster(:);
displaydataind_final_cell = displaydataind_final_cell(~cellfun(@isempty,displaydataind_final_cell));
displaydataind_final = zeros(sum(cellfun(@length,displaydataind_final_cell)),1);
counter = 1;
for ii = 1:length(displaydataind_final_cell)
    tempdisplaydataind_final_cell = displaydataind_final_cell{ii};
    displaydataind_final(counter:counter + length(tempdisplaydataind_final_cell) - 1) = ...
        tempdisplaydataind_final_cell(:);
    counter = counter + length(tempdisplaydataind_final_cell);
end
result = result(displaydataind_final);
save(fullfile(p,[f,'_O_bestisomers',e]),'result','scoreintdata');
end

function result_rev = recalcrevresult(result,displaydataind,allfragments,searchind,spectra,scannum,...
    denoisingoptCID,denoisingoptHCD,denoisingoptETD,denoisingoptETciD,denoisingoptEThcD,...
    scoreintdata,agpoptions)
result_rev = [];
ms2tol = scoreintdata.scoreoptions.ms2tol;
ms2tolunit = scoreintdata.scoreoptions.ms2tolunit;
resultsearchind = displaydataind(searchind,:);
allpepdata_rev = scoreintdata.prot;
for ii = size(resultsearchind,1):-1:1
    tempspectra = cell(1,size(resultsearchind,2));
    for jj = 1:size(resultsearchind,2)
        tempresult = result(resultsearchind(ii,jj));
        tempspectrum = spectra{scannum == tempresult.Scan};
        switch upper(result(resultsearchind(ii,jj)).Fragmode)
            case 'CID'
                methodind = strcmpi(scoreintdata.scoreoptions.analyzefragmode,'CID');
                tempspectrum = spectradenoising(tempspectrum,{'CID'},tempresult.Charge,...
                    tempresult.Theo,...
                    ms2tol(methodind),ms2tolunit{methodind},denoisingoptCID);
            case 'HCD'
                methodind = strcmpi(scoreintdata.scoreoptions.analyzefragmode,'HCD');
                tempspectrum = spectradenoising(tempspectrum,{'HCD'},tempresult.Charge,...
                    tempresult.Theo,...
                    ms2tol(methodind),ms2tolunit{methodind},denoisingoptHCD);
            case 'ETD'
                methodind = strcmpi(scoreintdata.scoreoptions.analyzefragmode,'ETD');
                tempspectrum = spectradenoising(tempspectrum,{'ETD'},tempresult.Charge,...
                    tempresult.Theo,...
                    ms2tol(methodind),ms2tolunit{methodind},denoisingoptETD);
            case 'ETCID'
                methodind = strcmpi(scoreintdata.scoreoptions.analyzefragmode,'ETCID');
                tempspectrum = spectradenoising(tempspectrum,{'ETCID'},tempresult.Charge,...
                    tempresult.Theo,...
                    ms2tol(methodind),ms2tolunit{methodind},denoisingoptETciD);
            case 'ETHCD'
                methodind = strcmpi(scoreintdata.scoreoptions.analyzefragmode,'ETHCD');
                tempspectrum = spectradenoising(tempspectrum,{'ETHCD'},tempresult.Charge,...
                    tempresult.Theo,...
                    ms2tol(methodind),ms2tolunit{methodind},denoisingoptEThcD);
        end
        tempspectra{jj} = tempspectrum;
    end
    revpepcomp = str2num(tempresult.ProteinID);
    revpepcomp(3:end) = flip(revpepcomp(3:end));
    revpepseq = scoreintdata.prot{revpepcomp(1)}{revpepcomp(2)};
    [p,g,m] = breakGlyPep(revpepseq);
    p.pep = flip(p.pep);
    if ~isempty(g)
        for jj = 1:length(g)
            g(jj).pos = length(p.pep) + 1 - g(jj).pos;
        end
    end
    if ~isempty(m)
        for jj = 1:length(m)
            m(jj).pos = length(p.pep) + 1 - m(jj).pos;
        end
    end
    revpepseq = joinGlyPep(p,g,m);
    allpepdata_rev{revpepcomp(1)}{revpepcomp(2)} = revpepseq;
    tempresult_rev = analyzegp({revpepcomp},result(ii).Theo,...
        [result(resultsearchind(ii,:)).Scan],[result(resultsearchind(ii,:)).Charge],tempspectra,...
        [result(resultsearchind(ii,:)).Retime],[result(resultsearchind(ii,:)).Expt],[result(resultsearchind(ii,:)).Mono],...
        [result(resultsearchind(ii,:)).Quant],1:size(displaydataind,2),1,allpepdata_rev,...
        scoreintdata.protname,1:length(allpepdata_rev),allfragments,...
        agpoptions);
    result_rev = [result_rev,tempresult_rev];
end
end