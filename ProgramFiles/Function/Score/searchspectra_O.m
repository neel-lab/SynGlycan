function searchspectra_O(input,msdata,triggerdata,allpepdata,statusreporthandles,...
    todo,allfragments,scoreoptions)

% This file sets up O-glycan data serach performed in the following
% 'searchogly.m' file. It:
% i) msdata: reads msdata file
% ii) agpoptions: sets up the scoring options
% iii) read PTM parameters from allfragments and allpepdata
% iv) sets up serach in doparallel and in single processor modes. This
% search is conducted by 'searchogly.m'
% v) The emerging 'result' are written into by 'savescoreresults.m' into input.outputfname


if ~isempty(statusreporthandles)
    statusstr = get(statusreporthandles.edit_statusreport, 'String');
end

% Read MSDATA and calculate precmass (precursor mass)
scannum = msdata.scannum;
charge = msdata.charge;
precmz = msdata.precursormz;
mslvl = msdata.mslvl;
fragmode = msdata.fragmode;
spectra = msdata.spectra;
retime = msdata.retime;
quant = msdata.parsedAUC;
allprecmz = msdata.allprecursormz;
exptmass = zeros(size(scannum));
for ii = 1:length(allprecmz)
    if ~isempty(allprecmz{ii})
        exptmass(ii) = (allprecmz{ii}(2) - 1.007825032) * charge(ii);
    end
end
precmass   = (precmz - 1.007825032).*charge;

% Get scoring options in the form agpoptions, for trigger mode and for
% regular mode
agpoptions.ms1tolunit = scoreoptions.ms1tolunit;
agpoptions.selectpeak = scoreoptions.selectpeak;
agpoptions.maxlag = scoreoptions.maxlag;
agpoptions.theofragopt_maxstublen = scoreoptions.maxstublen;
agpoptions.doparallel = scoreoptions.doparallel;
agpoptions.proceed_override = GlycoPATConstant.O_Proceed_Override;
if input.dohcdtrigger  % HCD trigger - supplemental activation experiments
    SASSO = triggerdata.SASSO;
    colnames = triggerdata.colnames;
    colselected = ismember(upper(colnames),upper(scoreoptions.analyzefragmode));
    colnames = colnames(colselected);
    SASSO = [SASSO(:,colselected),SASSO(:,end-1:end)];
    [~,ind1] = sort(colnames);
    [~,ind2] = sort(scoreoptions.analyzefragmode);
    [~,ind3] = sort(ind1);
    ind = ind2(ind3);
    SASSO(SASSO(:,2) == 0,:) = [];  % HCD THAT TRIGGERED NOTHING IS IGNORED!
    agpoptions.ms2tol = scoreoptions.ms2tol(ind);
    agpoptions.ms2tolunit = scoreoptions.ms2tolunit(ind);
    agpoptions.fragnum = scoreoptions.fragnum(ind,:);
    agpoptions.userpepiontyp = scoreoptions.userpepiontyp(ind);
    agpoptions.fragmode = colnames;
    agpoptions.theofragopt_monosacislabile = scoreoptions.monosacislabile(ind);
    agpoptions.theofragopt_simultaneousfrag = scoreoptions.simultaneousfrag(ind);
    agpoptions.theofragopt_addoxoniumion = scoreoptions.addoxoniumion(ind);
    for ii = 1:size(SASSO,1)
        for jj = 1:length(colnames)
            tempscanind = scannum == SASSO(ii,jj);
            fragmethodind = ismember(upper(scoreoptions.analyzefragmode),upper(colnames{jj}));
            tempms2tol = scoreoptions.ms2tol(fragmethodind);
            tempms2tolunit = scoreoptions.ms2tolunit{fragmethodind};
            spectra{tempscanind} = spectradenoising(spectra{tempscanind},colnames{jj},charge(tempscanind),...
                precmass(tempscanind),tempms2tol,tempms2tolunit,input);
        end
    end
else
    ms1searchscannum = scannum(mslvl == 2 &...
        ismember(fragmode,scoreoptions.analyzefragmode));
    for ii = 1:length(ms1searchscannum)
        tempscanind = scannum == ms1searchscannum(ii);
        tempfragmode = fragmode{tempscanind};
        tempms2tol = input.ms2tol(strcmpi(input.fragmode,tempfragmode));
        tempms2tolunit = input.ms2tolunit{strcmpi(input.fragmode,tempfragmode)};
        spectra{tempscanind} = spectradenoising(spectra{tempscanind},tempfragmode,charge(tempscanind),...
            precmass(tempscanind),tempms2tol,tempms2tolunit,input);
    end
    colnames = {'HCD'};
    agpoptions.ms2tol = scoreoptions.ms2tol;
    agpoptions.ms2tolunit = scoreoptions.ms2tolunit;
    agpoptions.fragnum = scoreoptions.fragnum;
    agpoptions.userpepiontyp = scoreoptions.userpepiontyp;
    agpoptions.fragmode = colnames;
    agpoptions.theofragopt_monosacislabile = scoreoptions.monosacislabile;
    agpoptions.theofragopt_simultaneousfrag = scoreoptions.simultaneousfrag;
    agpoptions.theofragopt_addoxoniumion = scoreoptions.addoxoniumion;
    triggerdata = [];
end

% Organize 'allfragments' data into PTM parameters
outputdir = input.outputdir;
outputfilename = input.outputfname;
[~,outputf,outputext] = fileparts(input.outputfname);
ptmfragnums = scoreoptions.fragnum(:,2);
ptmsortind = zeros(1,length(allpepdata.ptminfo.mod));
for ii = 1:length(allpepdata.ptminfo.mod)
    ptmsortind(ii) = find(ismember(allfragments.ptmseq,allpepdata.ptminfo.mod{ii}));
end
PTM.ptmseq = allfragments.ptmseq(ptmsortind);
PTM.ptmmass = allfragments.ptmmass(ptmsortind);
PTM.ptmmassfix = allfragments.ptmmassfix(ptmsortind);
PTM.ptmtype = allfragments.ptmtype(ptmsortind);
PTM.ptmfragsto = allfragments.ptmfragsto(ptmfragnums,ptmsortind);
PTM.ptmfragstublen = allfragments.ptmfragstublen(ptmfragnums,ptmsortind);
PTM.ptmuniind = allpepdata.ptminfo.uniind;  % PTM locator
PTM.ptmisglycan = allfragments.ptmisglycan(ptmsortind);
PTM.ptmisvar = allfragments.ptmisvar(ptmsortind);
[PTM.ptmuniind,ind] = sort(PTM.ptmuniind);
PTM.ptmseq = PTM.ptmseq(ind);
PTM.ptmmass = PTM.ptmmass(ind);
PTM.ptmmassfix = PTM.ptmmassfix(ind);
PTM.ptmtype = PTM.ptmtype(ind);
PTM.ptmfragsto  =PTM.ptmfragsto(:,ind);
PTM.ptmfragstublen = PTM.ptmfragstublen(:,ind);
PTM.ptmsize = cellfun(@(x) length(strfind(x,'{')),PTM.ptmseq);
PTM.ptmisglycan = PTM.ptmisglycan(ind);
PTM.ptmisvar = PTM.ptmisvar(ind);
PTM.isomers = PTM.ptmseq;

% setup for parallel computing
numprotperbatch = ceil(input.doparacomp * scoreoptions.numprotperworker);
totalprotnum = length(allpepdata.glypep);
[glycancombi,glycancombimass] = getglycancombi(PTM,scoreoptions.maxOglyonpep);

% Use 'fragOGpep.m' to organize pepfraginfo
for ii = 1:ceil(totalprotnum/numprotperbatch)
    protindse = [(ii-1)*numprotperbatch+1,min(ii*numprotperbatch,length(allpepdata.glypep))];
    [pepseq,pepcomp,pepmass,pepfrag_HCD,~,pepfrag_others,fragmode_others] = ...
        fragOGpep(allpepdata.glypep,protindse,agpoptions);  % pepcomp has been fixed to real values
    pepfraginfo.pepseq = pepseq;
    pepfraginfo.pepcomp = pepcomp;
    pepfraginfo.pepmass = pepmass;
    pepfraginfo.pepfrag_HCD = pepfrag_HCD;
    pepfraginfo.pepfrag_others = pepfrag_others;
    pepfraginfo.fragmode_others = fragmode_others;
%    save(fullfile(input.outputdir,[outputf,'_O_pepfrag',outputext]),'pepfraginfo');
    clear pepfraginfo;
    FASTAheadall = allpepdata.FASTAhead;

    % doparallel for trigger and HCD only runs. This basically calls 'searchogly.m'
    if scoreoptions.doparallel
        if input.dohcdtrigger
            SASSOind = 1:size(SASSO,1);
            dist_SASSOind = distributed(SASSOind);
            spmd
                [result_dist,newprotbatch_dist] = searchogly(pepseq,pepcomp,pepmass,FASTAheadall,...
                    pepfrag_HCD,pepfrag_others,fragmode_others,...
                    PTM,glycancombi,glycancombimass,scoreoptions,SASSO,colnames,...
                    agpoptions,scannum,spectra,charge,precmz,retime,exptmass,precmass,quant,...
                    getLocalPart(dist_SASSOind));
            end
        else  % HCD only
            HCDscannum = scannum(strcmpi(fragmode,'HCD'));
            HCDscannumind = 1:length(HCDscannum);
            dist_HCDscannumind = distributed(HCDscannumind);
            [result_dist,newprotbatch_dist] = searchogly(pepseq,pepcomp,pepmass,FASTAheadall,...
                pepfrag_HCD,pepfrag_others,fragmode_others,...
                PTM,glycancombi,glycancombimass,scoreoptions,HCDscannum(:),{'HCD'},...
                agpoptions,scannum,spectra,charge,precmz,retime,exptmass,precmass,quant,...
                getLocalPart(dist_HCDscannumind));
        end
        result = [];
        newprotbatch = {};
        oriprotid = [];
        for jj = 1:length(result_dist)
            tempresult_dist = result_dist{jj};
            result = [result,tempresult_dist];
            temporiprotid = zeros(length(tempresult_dist),2);
            for kk = 1:length(tempresult_dist)
                tempprotid = str2num(tempresult_dist(kk).ProteinID);
                temporiprotid(kk,:) = tempprotid(1:2);
            end
            oriprotid = [oriprotid;temporiprotid];
            newprotbatch = [newprotbatch;newprotbatch_dist{jj}];
        end
    else  % single thread for the same as above for trigger and HCD alone modes
        if input.dohcdtrigger
            SASSOind = 1:size(SASSO,1);
            [result,newprotbatch] = searchogly(pepseq,pepcomp,pepmass,FASTAheadall,...
                pepfrag_HCD,pepfrag_others,fragmode_others,...
                PTM,glycancombi,glycancombimass,scoreoptions,SASSO,colnames,...
                agpoptions,scannum,spectra,charge,precmz,retime,exptmass,precmass,quant,...
                SASSOind);
        else
            HCDscannum = scannum(strcmpi(fragmode,'HCD'));
            [result,newprotbatch] = searchogly(pepseq,pepcomp,pepmass,FASTAheadall,...
                pepfrag_HCD,pepfrag_others,fragmode_others,...
                PTM,glycancombi,glycancombimass,scoreoptions,HCDscannum(:),{'HCD'},...
                agpoptions,scannum,spectra,charge,precmz,retime,exptmass,precmass,quant,...
                1:length(HCDscannum));
        end
        oriprotid = zeros(length(result),2);
        for jj = 1:length(result)
            tempprotid = str2num(result(jj).ProteinID);
            oriprotid(jj,:) = tempprotid(1:2);
        end
    end

    % organize 'result' and write to file using 'savescoreresults.m'.
    % Output file name comes from 'ScoreAllSpectra.m'
    newglypep = cell(size(allpepdata.glypep));
    [uniprotind,~,uniprotindind] = unique(oriprotid(:,1));
    newprotid = oriprotid;
    for jj = 1:max(uniprotindind)
        tempprotbatch = newprotbatch(uniprotindind == jj);
        [unitempprotbatch,~,uniprotindsto] = unique(tempprotbatch,'stable');
        newglypep{uniprotind(jj)} = unitempprotbatch;
        newprotid(uniprotindind == jj,2) = uniprotindsto;
    end
    for jj = 1:length(result)
        tempprotid = str2num(result(jj).ProteinID);
        tempprotid(2) = newprotid(jj,2);
        result(jj).ProteinID = num2str(tempprotid);
    end
    allpepdata.glypep = newglypep;
    savescoreresults(result,scoreoptions,input,PTM,allpepdata,triggerdata,...
        outputdir,outputfilename,protindse,'save');
    
    if ~isempty(statusreporthandles)
        analytime = toc;
        statusstr = [['Group ',num2str(ii),': Protein No.',num2str(protindse(1)),...
            ' to ',num2str(protindse(2)),' finished, time = ',num2str(analytime),' sec.'];...
            [num2str(totalprotnum-protindse(2)),' remaining.'];'Saving results...';statusstr];
        set(statusreporthandles.edit_statusreport, 'String', statusstr);
        pause(0.0001);
    end
end

if ~isempty(statusreporthandles)
    statusstr = ['Running isomeric glycopeptide calculations...';statusstr];
    set(statusreporthandles.edit_statusreport, 'String', statusstr);
    pause(0.0001);
end

% findbestglypep_O(input,PTM,agpoptions,scannum,spectra,allpepdata); Getting
% rid of findbestisomers_O

% findbestglypep_O(input,allfragments,agpoptions,scannum,spectra,allpepdata);

if ~isempty(statusreporthandles)
    statusstr = ['All calculations finished';statusstr];
    set(statusreporthandles.edit_statusreport, 'String', statusstr);
    pause(0.0001);
end
end