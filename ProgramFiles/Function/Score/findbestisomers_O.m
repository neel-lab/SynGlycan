function findbestisomers_O(scorefilename,ispara)
if isstruct(scorefilename)
    scorefilename = fullfile(scorefilename.folder,scorefilename.name);
elseif ischar(scorefilename)
    % INTENTIONALLY LEFT BLANK
end
load(scorefilename,'result','scoreintdata');
load([scorefilename(1:end-7),'fragment.mat'],'allfragments');
scoreintdata_withfrag = scoreintdata;
scoreintdata_withfrag.ptmfragsto = allfragments.ptmfragsto;
scoreintdata_withfrag.ptmfragstublen = allfragments.ptmfragstublen;
scoreintdata_withfrag.ptmmass = allfragments.ptmmass;
scoreintdata_withfrag.ptmtype = allfragments.ptmtype;
scoreintdata_withfrag.ptmseq = allfragments.ptmseq;
if exist(scoreintdata.sliminput.exptdata,'file')
    load(scoreintdata.sliminput.exptdata,'msdata_QA');
else
    [f,p] = uigetfile('*.mat','Load experiment data file');
    load(fullfile(p,f),'msdata_QA');
end
allspectra = msdata_QA.spectra;
allchg = msdata_QA.charge;
allscan = msdata_QA.scannum;
allprecmz = msdata_QA.precursormz;
allprecmass = (allprecmz - 1.007825032) .* allchg;
resultid = cell(size(result));
for ii = 1:length(result)
    resultid{ii} = [num2str(result(ii).Scan),'_',result(ii).ProteinID];
end
[~,uniind,~] = unique(resultid);
result = result(uniind);
protids = {result.ProteinID};
scans = [result.Scan];

displaydataind = grouptriggeredresult(scoreintdata.SASSO,protids,scans,'combineisomer');
originalresult = result;
clear result;
displaydataind_ind = 1:size(displaydataind,1);
if ispara
    dist_displaydataind_ind = distributed(displaydataind_ind);
    spmd
        result_dist = calcsgp_O(displaydataind,getLocalPart(dist_displaydataind_ind),originalresult,...
        allscan,allchg,allprecmass,allspectra,scoreintdata_withfrag);
    end
    result = [];
    for ii = 1:length(result_dist)
        result = [result,result_dist{ii}];
    end
else
    result = calcsgp_O(displaydataind,displaydataind_ind,originalresult,...
        allscan,allchg,allprecmass,allspectra,scoreintdata_withfrag);
end
save([scorefilename(1:end-4),'_bestisomers.mat'],'result','scoreintdata');
end