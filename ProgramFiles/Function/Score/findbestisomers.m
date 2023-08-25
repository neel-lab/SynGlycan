function findbestisomers(input,msdata_QA,allfragments,rankingmethod,allpepdata)
% input:
% doparacomp
% outputfname
% outputdir
% numprotperworker
% fbioptions:
% allpepdata
% nummindiagionfragmode




% wtbar = waitbar(0,'Calculating');
isogpEScutoff = GlycoPATConstant.EScutoff_fbi;
pepcovcutoff = GlycoPATConstant.Pepcovcutoff_fbi;

numprotperbatch = input.doparacomp * input.numprotperworker;
% allpepdata = fbioptions.allpepdata;
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
    tempoutputfilename = fullfile(p,[f,'_',num2str(protindse(1)),'_',num2str(protindse(2)),e]);
    if isfield(input,'peporglycopep')
        tempoutputfilename = fullfile(p,[f,'_N_glycopep','_',num2str(protindse(1)),'_',num2str(protindse(2)),e]);
    end
    if ismember(tempoutputfilename,existingdatafilenames)
        if ii == 1
            load(tempoutputfilename,'result','scoreintdata');
        else
            load(tempoutputfilename,'result');
        end
        accumulatedresult = [accumulatedresult,result];
    else
        tempoutputfilename_regexppattern = [f,'_',num2str(protindse(1)),'_',...
            num2str(protindse(2)),'_Part_[0-9]+of[0-9]+',e];
        fileexists = ~cellfun(@isempty,regexp(existingdatafilenames,tempoutputfilename_regexppattern));
        temppresearchfilenames = existingdatafilenames(fileexists);
        for jj = 1:length(temppresearchfilenames)
            tempoutputfilename = temppresearchfilenames{jj};
            if (ii == 1)  && (jj == 1)
                load(tempoutputfilename,'result','scoreintdata');
            else
                load(tempoutputfilename,'result');
            end
            accumulatedresult = [accumulatedresult,result];
        end
    end
end
result = accumulatedresult;
clear accumulatedresult
allspectra = msdata_QA.spectra;
allchg = msdata_QA.charge;
allscan = msdata_QA.scannum;
resultid = cell(size(result));
for ii = 1:length(result)
    resultid{ii} = [num2str(result(ii).Scan),'_',result(ii).ProteinID];
end
[~,uniind,~] = unique(resultid);
result = result(uniind);
protids = {result.ProteinID};
scans = [result.Scan];
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
if scoreintdata.sliminput.dohcdtrigger
    displaydataind = grouptriggeredresult(scoreintdata.SASSO,protids,scans,'regular'); %172 groups of scans
else
    displaydataind =  (1:length(result))';
end
displayoptions.scoreintdata = scoreintdata;
displayoptions.combineisomer = false;
if scoreintdata.scoreoptions.isHCDtrigger
    [CombiES,~] = calculateCombiES(result,displaydataind,...
        scoreintdata.colnames);
    switch rankingmethod
        case 1
            escutoff = fdrfilter(result,displayoptions,...
                displaydataind,.01,'CombiES');
            rowkeepind = CombiES >= escutoff;
            for ii = 1:size(displaydataind,1)
                tempenscore = [result(displaydataind(ii,:)).Enscore];
                if any(tempenscore == 0)
                    rowkeepind(ii) = false;
                end
            end
        case 2
            ScanHCD = [result(displaydataind(:,1)).Scan];
            [~,~,ind2] = unique(ScanHCD);
            rowkeepind = false(size(displaydataind,1),1);
            for ii = 1:max(ind2)
                indices = find(ind2 == ii);
                tempCombiES = CombiES(indices);
                topCombiES = sort(tempCombiES,'descend');
                if length(topCombiES) > 1
                    if topCombiES(1) - topCombiES(2) <= 0.05
                        rowkeepind(indices((tempCombiES >= isogpEScutoff) & (tempCombiES >= topCombiES(2)))) = true;
                    else
                        rowkeepind(indices((tempCombiES >= isogpEScutoff) & (tempCombiES == max(tempCombiES)))) = true;
                    end
                else
                    rowkeepind(indices((tempCombiES >= isogpEScutoff) & (tempCombiES == max(tempCombiES)))) = true;
                end
            end
            rowkeepind_pep = false(size(displaydataind,1),1);
            for ii = 1:size(displaydataind,1)
                rowkeepind_pep(ii) = max([result(displaydataind(ii,:)).Pepcov]) >= pepcovcutoff;
            end
            rowkeepind = rowkeepind & rowkeepind_pep;
    end
    displaydataind = displaydataind(rowkeepind,:);
    displaydataind = displaydataind';
    result = result(displaydataind(:));
    displaydataind = reshape(1:length(result),length(scoreintdata.colnames),...
        length(result)/length(scoreintdata.colnames));
    displaydataind = displaydataind';
else
    escutoff = fdrfilter(result,displayoptions,...
        displaydataind,.01,'Every hit');
    rowkeepind = [result.Enscore] >= escutoff;
    displaydataind = displaydataind(rowkeepind,:);
    displaydataind = displaydataind';
    result = result(displaydataind(:));
    %%%%%%%%%%%%%%%% deleting reshape
    %displaydataind = reshape(1:length(result),length(scoreintdata.colnames),...
        %length(result)/length(scoreintdata.colnames));
    displaydataind = displaydataind';
end

if isfield(input,'peporglycopep')
    save(fullfile(p,[f,'_N_glycopep_resultbeforecalc.mat']),'result');  
else
%    save(fullfile(p,[f,'_resultbeforecalc.mat']),'result');
end

sgps = {result.SGP};
scans = [result.Scan];
allptmseq = allfragments.ptmseq;
ptmisomers = scoreintdata.ptmisomers;
allfragment_ptmisomerind = zeros(size(ptmisomers));
for ii = 1:size(ptmisomers,1)
    for jj = 1:size(ptmisomers,2)
        if ~isempty(ptmisomers{ii,jj})
            allfragment_ptmisomerind(ii,jj) = find(ismember(allptmseq,ptmisomers{ii,jj}));
        end
    end
end
allfragmentkeepind = false(size(allfragments.ptmseq));
for ii = 1:length(result)
    [~,g,m] = breakGlyPep(result(ii).SGP);
    if ~isempty(g)
        for jj = 1:length(g)
            tempptm = g(jj).struct;
            ptmisomerlocation = strcmpi(ptmisomers(:,1),tempptm);
            tempallfragmentkeepind = allfragment_ptmisomerind(ptmisomerlocation,:);
            tempallfragmentkeepind = tempallfragmentkeepind(tempallfragmentkeepind > 0);
            allfragmentkeepind(tempallfragmentkeepind) = true;
        end
    end
    if ~isempty(m)
        for jj = 1:length(m)
            tempptm = m(jj).struct;
            ptmisomerlocation = strcmpi(ptmisomers(:,1),tempptm);
            tempallfragmentkeepind = allfragment_ptmisomerind(ptmisomerlocation,:);
            tempallfragmentkeepind = tempallfragmentkeepind(tempallfragmentkeepind > 0);
            allfragmentkeepind(tempallfragmentkeepind) = true;
        end
    end
end
for ii = 1:length(allptmseq)
    if ~allfragmentkeepind(ii)
        for jj = 1:size(allfragments.ptmfragsto,1)
            allfragments.ptmfragsto{jj,ii} = {};
            allfragments.ptmfragstublen{jj,ii} = {};
        end
    end
end

numpw = scoreintdata.sliminput.doparacomp;

doparallel = numpw > 1;
searchind = 1:size(displaydataind,1);
allspectra = allspectra(ismember(allscan,[result.Scan]));
allchg = allchg(ismember(allscan,[result.Scan]));
allscan = intersect(allscan,[result.Scan]);
clear msdata
if scoreoptions.isHCDtrigger
    if ~isfield(scoreintdata,'colnames')
        scoreintdata.colnames = scoreintdata.scoreoptions.analyzefragmode;
    end
    skip=1;
    if skip==0
%    if doparallel
        dist_searchind = distributed(searchind);
        isomericglypepresult = [];
        spmd
            dist_isomericglypepresult = calcsgp(result,sgps,scans,displaydataind,allspectra,allscan,allchg,...
                denoisingoptHCD,denoisingoptCID,denoisingoptETD,denoisingoptETciD,denoisingoptEThcD,scoreintdata,allfragments,...
                getLocalPart(dist_searchind));
        end
        for ii = 1:length(dist_isomericglypepresult)
            isomericglypepresult = [isomericglypepresult,dist_isomericglypepresult{ii}];
        end
    else
        isomericglypepresult = calcsgp(result,sgps,scans,displaydataind,allspectra,allscan,allchg,...
            denoisingoptHCD,denoisingoptCID,denoisingoptETD,denoisingoptETciD,denoisingoptEThcD,scoreintdata,allfragments,searchind);
    end
else
%{
    if doparallel
        dist_searchind = distributed(searchind);
        isomericglypepresult = [];
        spmd
            dist_isomericglypepresult = calcsgp_legacy(result,sgps,scans,displaydataind,allspectra,allscan,allchg,...
                denoisingoptHCD,denoisingoptCID,denoisingoptEThcD,scoreintdata,allfragments,...
                getLocalPart(dist_searchind));
        end
        for ii = 1:length(dist_isomericglypepresult)
            isomericglypepresult = [isomericglypepresult,dist_isomericglypepresult{ii}];
        end
    else
%}
    isomericglypepresult = calcsgp_legacy(result,sgps,scans,displaydataind,allspectra,allscan,allchg,...
        denoisingoptHCD,denoisingoptCID,denoisingoptEThcD,scoreintdata,allfragments,searchind);
%   end
end
result = isomericglypepresult;
if isfield(input,'peporglycopep')
    save(fullfile(p,[f,'_N_glycopep_bestisomers.mat']),'result','scoreintdata');
else
    save(fullfile(p,[f,'_bestisomers.mat']),'result','scoreintdata');
end
        
end