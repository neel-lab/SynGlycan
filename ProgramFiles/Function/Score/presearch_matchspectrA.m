function result = presearch_matchspectrA(msdata,matched_gpcomp,matched_gpmass,...
    matched_scanind,protbatch,FASTAbatch,triggerdata,protind,PTM,scoreoptions)
% PRESEARCH_MATCHSPECTRA: match experiment spectra against glycopeptides.
%
% Syntax:
% result = presearch_matchspectrA(msdata,matched_gpcomp,matched_gpmass,...
%     matched_scanind,protbatch,FASTAbatch,triggerdata,protind,PTM,scoreoptions)
%
% Input:
% msdata: structure? Experiment data.
% matched_gpcomp: n x 1 cell array of strings. Composition of matched
%     candidates.
% matched_gpmass: n x 1 double. Theoretical mass of matched candidates.
% matched_scanind: n x 1 double. Index of matched spectrum.
%     For above 3 inputs see FINDMS1MATCH for detail.
% protbatch: n x 1 cell array of m x 1 cell array of strings. Candidate
%     peptide backbones.
% FASTAbatch: n x 1 cell array of strings. Protein FASTA names.
% triggerdata: structure. Scan association for triggered experiments. See
%     SCOREALLSPECTRA for detail.
% protind: 1 x n double. Protein index.
% PTM: structure. PTMs of glycopeptides.
% scoreoptions: structure. Scoring options. See SCOREALLSPECTRA for detail.
%
% Output:
% result: 1 x n structure. Scoring results. See MATCH1BY1 for detail.
%
% Note:
% N/A
%
% Example:
% N/A
%
% Children function:
% N/A
%

result = [];
%% UNPACK DATA & OPTIONS
isHCDtrigger = scoreoptions.isHCDtrigger;
analyzefragmode = scoreoptions.analyzefragmode;
fragnum = scoreoptions.fragnum;
doparallel = scoreoptions.doparallel;
parallelactive = scoreoptions.parallelactive;

%% PREPARE EXPT DATA
scannum = msdata.scannum;
charge = msdata.charge;
precmz = msdata.precursormz;
fragmode = msdata.fragmode;
precmass = (precmz - 1.007825032).*charge;
exptmz = zeros(size(precmz));
for ii = 1:length(precmz)
    tempprecmz = msdata.allprecursormz{ii};
    if ~isempty(tempprecmz)
        exptmz(ii) = tempprecmz(2);
    end
end
exptmass = (exptmz - 1.007825032).*charge;

%% START ANALYSIS
if isHCDtrigger
    SASSO = triggerdata.SASSO;
    colnames = triggerdata.colnames;
    colselected = ismember(upper(colnames),upper(analyzefragmode));
    colnames = colnames(colselected);
    SASSO = [SASSO(:,colselected),SASSO(:,end-1:end)];
    [~,ind1] = sort(colnames);
    [~,ind2] = sort(analyzefragmode);
    [~,ind3] = sort(ind1);
    ind = ind2(ind3);
    denoisingoptions.fragmode = scoreoptions.analyzefragmode;
    denoisingoptions.cutoffmed = scoreoptions.cutoffmed;
    denoisingoptions.fracmax = scoreoptions.fracmax;
    denoisingoptions.minmaxmz = scoreoptions.minmaxmz;
    
    % See ANALYZEGP and SCOREALLSPECTRA for details about these options below.
    agpoptions.ms1tolunit = scoreoptions.ms1tolunit;
    agpoptions.ms2tol = scoreoptions.ms2tol(ind);
    agpoptions.ms2tolunit = scoreoptions.ms2tolunit(ind);
    agpoptions.fragnum = fragnum(ind,:);
    agpoptions.fragmode = colnames;
    agpoptions.selectpeak = scoreoptions.selectpeak;
    agpoptions.maxlag = scoreoptions.maxlag;
    agpoptions.theofragopt_monosacislabile = scoreoptions.monosacislabile(ind);
    agpoptions.theofragopt_simultaneousfrag = scoreoptions.simultaneousfrag(ind);
    agpoptions.theofragopt_addoxoniumion = scoreoptions.addoxoniumion(ind);
    agpoptions.theofragopt_maxstublen = scoreoptions.maxstublen;
    PTM.ptmfragsto = PTM.ptmfragsto(ind,:);
    PTM.ptmfragstublen = PTM.ptmfragstublen(ind,:);
    searchscannum = SASSO(matched_scanind,1:end-2);
    [unisearchscannum,~,unisearchscannumind] = unique(searchscannum);
    unisearchscannumind = reshape(unisearchscannumind,size(searchscannum,1),...
        size(searchscannum,2));
    extractind = ismember(scannum,unisearchscannum);
    spectrasto = msdata.spectra(extractind);
    chargesto = charge(extractind);
    retimesto = msdata.retime(extractind);
    exptmasssto = exptmass(extractind);
    monomasssto = precmass(extractind);
    quantsto = msdata.parsedAUC(extractind);
    fragmodesto = msdata.fragmode(extractind);
    ms2tols = zeros(size(spectrasto));
    ms2tolunits = cell(size(spectrasto));
    allfragmodes = agpoptions.fragmode;
    allms2tols = agpoptions.ms2tol;
    allms2tolunits = agpoptions.ms2tolunit;
    for ii = 1:length(spectrasto)
        fragmodeind = ismember(upper(allfragmodes),upper(fragmodesto{ii}));
        ms2tols(ii) = allms2tols(fragmodeind);
        ms2tolunits{ii} = allms2tolunits(fragmodeind);
    end
    spectrasto = spectradenoising(spectrasto,fragmodesto,chargesto,...
        monomasssto,ms2tols,ms2tolunits,denoisingoptions);
    clear matched_scanind searchscannum triggerdata SASSO ind* extractind ms2tols ms2tolunits
    clear msdata scannum charge precmz precmass exptmz exptmass fragmentsto
    unisearchscannumindind = 1:size(unisearchscannumind,1);
    if doparallel
        dist_matched_gpcomp = distributed(matched_gpcomp);
        dist_matched_gpmass = distributed(matched_gpmass);
        dist_unisearchscannumindind = distributed(unisearchscannumindind);
        spmd
            dist_tempresult = presearch_analyzegp(getLocalPart(dist_matched_gpcomp),...
                getLocalPart(dist_matched_gpmass),unisearchscannum,chargesto,spectrasto,...
                retimesto,exptmasssto,monomasssto,quantsto,unisearchscannumind,...
                getLocalPart(dist_unisearchscannumindind),protbatch,FASTAbatch,protind,PTM,...
                agpoptions);
        end
        for ii = 1:length(dist_tempresult)
            result = [result,dist_tempresult{ii}];
        end
        clear dist_*
    else
        result = presearch_analyzegp(matched_gpcomp,matched_gpmass,unisearchscannum,chargesto,spectrasto,...
            retimesto,exptmasssto,monomasssto,quantsto,unisearchscannumind,...
            unisearchscannumindind,protbatch,FASTAbatch,protind,PTM,...
            agpoptions);
    end
else  % regular
    agpoptions = scoreoptions;
    agpoptions.fragmode = scoreoptions.analyzefragmode;
    agpoptions.theofragopt_monosacislabile = scoreoptions.monosacislabile;
    agpoptions.theofragopt_simultaneousfrag = scoreoptions.simultaneousfrag;
    agpoptions.theofragopt_addoxoniumion = scoreoptions.addoxoniumion;
    agpoptions.theofragopt_maxstublen = scoreoptions.maxstublen;
    mslvl = msdata.mslvl;
    searchscannum_ori = scannum(mslvl == 2 &...
        ismember(upper(msdata.fragmode),upper(scoreoptions.analyzefragmode)));
    searchfragmode = fragmode(mslvl == 2 &...
        ismember(upper(msdata.fragmode),upper(scoreoptions.analyzefragmode)));
    matchedsearchscannum = searchscannum_ori(matched_scanind);
    matchedsearchfragmode = searchfragmode(matched_scanind);
    glypepcomplength = cellfun(@length,matched_gpcomp);
    matched_gpcomp_ori = matched_gpcomp;
    matched_gpmass_ori = matched_gpmass;
    matched_glypepcomp = zeros(length(matched_gpcomp_ori),max(glypepcomplength));
    for ii = 1:length(matched_gpcomp)
        matched_glypepcomp(ii,1:glypepcomplength(ii)) = matched_gpcomp{ii};
    end
    [~,matched_uniglypepind1,matched_uniglypepind2] = unique(matched_glypepcomp,'stable','rows');
    matched_gpcomp = {};
    matched_gpmass = [];
    searchscannum = [];
    for ii = 1:max(matched_uniglypepind2)
        tempmatchedsearchscannum = matchedsearchscannum(matched_uniglypepind2 == ii);
        tempmatchedsearchfragmode = matchedsearchfragmode(matched_uniglypepind2 == ii);
        tempmatched_gpcomp = matched_gpcomp_ori(matched_uniglypepind1(ii));
        tempmatched_gpmass = matched_gpmass_ori(matched_uniglypepind1(ii));
        numscanperfrag = zeros(1,length(analyzefragmode));
        for jj = 1:length(analyzefragmode)
            numscanperfrag(jj) = sum(strcmpi(tempmatchedsearchfragmode,analyzefragmode{jj}));
        end
        tempsearchscannum = zeros(max(numscanperfrag),length(analyzefragmode));
        for jj = 1:length(analyzefragmode)
            thistempsearchscannum = tempmatchedsearchscannum(strcmpi(tempmatchedsearchfragmode,...
                analyzefragmode{jj}));
            tempsearchscannum(1:length(thistempsearchscannum),jj) = thistempsearchscannum(:);
        end
        searchscannum = [searchscannum;tempsearchscannum];
        matched_gpcomp = [matched_gpcomp;repmat(tempmatched_gpcomp,max(numscanperfrag),1)];
        matched_gpmass = [matched_gpmass;repmat(tempmatched_gpmass,max(numscanperfrag),1)];
    end
    [unisearchscannum,~,unisearchscannumind] = unique(searchscannum);
    if ismember(0,unisearchscannum)
        unisearchscannum = unisearchscannum(2:end);
        unisearchscannumind = unisearchscannumind - 1;
    end
    unisearchscannumind = reshape(unisearchscannumind,size(searchscannum,1),...
        size(searchscannum,2));
    extractind = ismember(scannum,unisearchscannum);
    spectrasto = msdata.spectra(extractind);
    chargesto = charge(extractind);
    retimesto = msdata.retime(extractind);
    exptmasssto = exptmass(extractind);
    monomasssto = precmass(extractind);
    quantsto = msdata.parsedAUC(extractind);
    fragmodesto = msdata.fragmode(extractind);
    ms2tols = zeros(size(spectrasto));
    ms2tolunits = cell(size(spectrasto));
    allfragmodes = agpoptions.fragmode;
    allms2tols = agpoptions.ms2tol;
    allms2tolunits = agpoptions.ms2tolunit;
    for ii = 1:length(spectrasto)
        fragmodeind = ismember(upper(fragmodesto{ii}),upper(allfragmodes));
        ms2tols(ii) = allms2tols(fragmodeind);
        ms2tolunits{ii} = allms2tolunits(fragmodeind);
    end
    denoisingoptions.fragmode = scoreoptions.analyzefragmode;
    denoisingoptions.cutoffmed = scoreoptions.cutoffmed;
    denoisingoptions.fracmax = scoreoptions.fracmax;
    denoisingoptions.minmaxmz = scoreoptions.minmaxmz;
    spectrasto = spectradenoising(spectrasto,fragmodesto,chargesto,...
        monomasssto,ms2tols,ms2tolunits,denoisingoptions);
    clear matched_scanind searchscannum triggerdata SASSO ind* extractind ms2tols ms2tolunits
    clear msdata scannum charge precmz precmass exptmz exptmass fragmentsto
    unisearchscannumindind = 1:size(unisearchscannumind,1);
    if doparallel
        dist_matched_gpcomp = distributed(matched_gpcomp);
        dist_matched_gpmass = distributed(matched_gpmass);
        dist_unisearchscannumindind = distributed(unisearchscannumindind);
        spmd
            dist_tempresult = presearch_analyzegp(getLocalPart(dist_matched_gpcomp),...
                getLocalPart(dist_matched_gpmass),unisearchscannum,chargesto,spectrasto,...
                retimesto,exptmasssto,monomasssto,quantsto,unisearchscannumind,...
                getLocalPart(dist_unisearchscannumindind),protbatch,FASTAbatch,protind,PTM,...
                agpoptions);
        end
        for ii = 1:length(dist_tempresult)
            result = [result,dist_tempresult{ii}];
        end
        clear dist_*
    else
        result = presearch_analyzegp(matched_gpcomp,matched_gpmass,unisearchscannum,chargesto,spectrasto,...
            retimesto,exptmasssto,monomasssto,quantsto,unisearchscannumind,...
            unisearchscannumindind,protbatch,FASTAbatch,protind,PTM,...
            agpoptions);
    end
end
end