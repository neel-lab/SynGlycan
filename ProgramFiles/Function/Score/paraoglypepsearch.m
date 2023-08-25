function paraoglypepsearch(msdata,input)
pepfile = input.pepfile;
exptdata = input.exptdata;
scoreoptions.analyzefragmode = input.fragmode;
scoreoptions.ms1tol = input.ms1tol;
scoreoptions.ms1tolunit = input.ms1tolunit;
scoreoptions.ms2tol = input.ms2tol;
scoreoptions.ms2tolunit = input.ms2tolunit;
scoreoptions.minmaxmz = input.minmaxmz;
scoreoptions.maxlag = input.maxlag;
scoreoptions.cutoffmed = input.cutoffmed;
scoreoptions.fracmax = input.fracmax;
scoreoptions.fragnum = input.fragnum;
scoreoptions.selectpeak = input.selectpeak;
scoreoptions.assemblegpmode = 1;
scoreoptions.numprotperworker = input.numprotperworker;
scoreoptions.resumeprevsearch = input.resumeprevsearch;
scoreoptions.excludenglysequon = input.excludenglysequon;

% Special options
scoreoptions.doparallel = false;
scoreoptions.monosacislabile = input.monosacislabile;
scoreoptions.simultaneousfrag = input.simultaneousfrag;
scoreoptions.addoxoniumion = input.addoxoniumion;
scoreoptions.isHCDtrigger = input.dohcdtrigger;
scoreoptions.maxstublen = input.maxstublen;
scoreoptions.presearch = input.presearch;
scoreoptions.maxglycannum = input.maxglycannum;
scoreoptions.keeptopn = input.keeptopnglypep;

agpoptions = scoreoptions;

allfragments = input.allfragments;
defaultfragiontyp = 'bcyzi';
fragiontyp = input.fragiontyp;
userpepiontyp = cell(size(fragiontyp,1),1);
for ii = 1:length(userpepiontyp)
    userpepiontyp{ii} = defaultfragiontyp(logical(fragiontyp(ii,:)));
end
scoreoptions.userpepiontyp = userpepiontyp;
allpepdata = digestfileanaly(pepfile,1);
glypep = allpepdata.glypep;
digestable_list = false(size(glypep));
ptmtyp = allpepdata.ptminfo.uniind;
ptmseq = allpepdata.ptminfo.mod;
[unityp,ind] = unique(ptmtyp);
isglycan = false(size(unityp));
for ii = 1:length(ind)
    isglycan(ii) = any(strfind(ptmseq{ind(ii)},'{'));  % Glycans are recognized by "{"
end
glycanind = unityp(isglycan);
for ii = 1:length(glypep)
    thisprot = glypep{ii};
    keepglypep = false(size(thisprot));
    % Identify the PTM marker on each glypep, if peptide does
    %     not carry glycan, abandon this one (Case 1: glypep only)
    for jj = 1:length(keepglypep)
        glypepseq = thisprot{jj};
        ptmtypes = regexp(glypepseq,'{[0-9]+}','match');
        ptmind = zeros(size(ptmtypes));  % Find the markers, check if they are glycans.
        for kk = 1:length(ptmtypes)
            ptmind(kk) = str2double(ptmtypes{kk}(2:end-1));
            if ismember(ptmind(kk),glycanind)
                keepglypep(jj) = true;
            end
        end
        oglycosylationpos = sort([strfind(glypepseq,'S'),strfind(glypepseq,'T')],'descend');
        tempglypepseq = glypepseq;
        for kk = 1:length(oglycosylationpos)
            thisoglypos = oglycosylationpos(kk);
            if (thisoglypos + 1 <= length(glypepseq)) && (~strcmpi(tempglypepseq(thisoglypos + 1),'{'))
                tempglypepseq = [tempglypepseq(1:thisoglypos),'{2}',tempglypepseq(thisoglypos+1:end)];
            end
        end
        thisprot{jj} = tempglypepseq;
        [p,~,~] = breakGlyPep(glypepseq);
        pepseq = p.pep;
        if any(regexp(pepseq,'N[^P][ST]'))
            keepglypep(jj) = false;
        end
    end
    
    if any(keepglypep)
        digestable_list(ii) = true;  % If a protein contains no glypep it's rejected
        glypep{ii} = unique(thisprot(keepglypep));
    end
end
allpepdata.glypep = glypep(digestable_list);
FASTAhead = allpepdata.FASTAhead;
allpepdata.FASTAhead = FASTAhead(digestable_list);
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

[path,filename,~] = fileparts(exptdata);
agpoptions.ms1tolunit = scoreoptions.ms1tolunit;
agpoptions.selectpeak = scoreoptions.selectpeak;
agpoptions.maxlag = scoreoptions.maxlag;
agpoptions.theofragopt_maxstublen = scoreoptions.maxstublen;

if input.dohcdtrigger  % HCD trigger - supplemental activation experiments
    sassofilename = fullfile(path,[filename,'_trigger.mat']);
    triggerdata = load(sassofilename);
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
    agpoptions.fragmode = colnames;
    agpoptions.theofragopt_monosacislabile = scoreoptions.monosacislabile;
    agpoptions.theofragopt_simultaneousfrag = scoreoptions.simultaneousfrag;
    agpoptions.theofragopt_addoxoniumion = scoreoptions.addoxoniumion;
    triggerdata = [];
end

outputdir = input.outputdir;
outputfilename = input.outputfname;
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

allfragments.ptmseq = PTM.ptmseq;
allfragments.ptmmass = PTM.ptmmass;
allfragments.ptmmassfix = PTM.ptmmassfix;
allfragments.ptmtype = PTM.ptmtype;
allfragments.ptmfragsto = PTM.ptmfragsto;
allfragments.ptmfragstublen = PTM.ptmfragstublen;
allfragments.ptmisglycan = PTM.ptmisglycan;
allfragments.ptmisvar = PTM.ptmisvar;
newfragfilename = fullfile(outputdir,[outputfilename(1:end-4),'_fragment.mat']);
save(newfragfilename,'allfragments');
input.allfragmentspath = newfragfilename;

numprotperbatch = ceil(input.doparacomp * scoreoptions.numprotperworker);
totalprotnum = length(allpepdata.glypep);
[glycancombi,glycancombimass] = getglycancombi(PTM,scoreoptions.maxglycannum);

for ii = 1:ceil(totalprotnum/numprotperbatch)
    protindse = [(ii-1)*numprotperbatch+1,min(ii*numprotperbatch,length(allpepdata.glypep))];
    protbatch = allpepdata.glypep(protindse(1):protindse(2));
    FASTAbatch = allpepdata.FASTAhead(protindse(1):protindse(2));
    protind = protindse(1):protindse(2);
    result = [];
    for jj = 1:length(protbatch)
        thisprotbatch = protbatch{jj};
        tempnewprotbatch = {};
        tempnewprotbatchindoffset = 0;
        thisprotbatchmass = cellfun(@pepMW,thisprotbatch);
        thisprotbatch_pepfrag = cell(size(thisprotbatch));
        if scoreoptions.doparallel
            parfor kk = 1:length(thisprotbatch)
                thisprotbatch_pepfrag{kk} = fragpep(thisprotbatch{kk},1,'bcyz','');
            end
            if input.dohcdtrigger
                SASSOind = 1:size(SASSO,1);
                dist_SASSOind = distributed(SASSOind);
                spmd
                    [result_dist,newprotbatch_dist] = spmd_searchogly(thisprotbatch,thisprotbatch_pepfrag,...
                        FASTAbatch{jj},thisprotbatchmass,PTM,glycancombi,glycancombimass,scoreoptions,SASSO,colnames,...
                        agpoptions,scannum,spectra,charge,precmz,retime,exptmass,precmass,quant,getLocalPart(dist_SASSOind));
                end
                for kk = 1:length(result_dist)
                    tempresult_dist = result_dist{kk};
                    tempnewprotbatch_dist = newprotbatch_dist{kk};
                    for ll = 1:length(tempresult_dist)
                        tempresult_dist(ll).ProteinID(1:2) = [protind(jj),...
                            tempresult_dist(ll).ProteinID(2) + tempnewprotbatchindoffset];
                    end
                    tempnewprotbatchindoffset = tempnewprotbatchindoffset + length(tempnewprotbatch_dist);
                    result = [result,tempresult_dist];
                    tempnewprotbatch = [tempnewprotbatch;tempnewprotbatch_dist];
                end
            else  % HCD only
                dist_SASSOind = distributed(1:length(ms1searchscannum));
                SASSO = ms1searchscannum(:);
                spmd
                    [result_dist,newprotbatch_dist] = spmd_searchogly(thisprotbatch,thisprotbatch_pepfrag,...
                        FASTAbatch{jj},thisprotbatchmass,PTM,glycancombi,glycancombimass,scoreoptions,SASSO,colnames,...
                        agpoptions,scannum,spectra,charge,precmz,retime,exptmass,precmass,quant,getLocalPart(dist_SASSOind));
                end
                for kk = 1:length(result_dist)
                    tempresult_dist = result_dist{kk};
                    tempnewprotbatch_dist = newprotbatch_dist{kk};
                    for ll = 1:length(tempresult_dist)
                        tempresult_dist(ll).ProteinID(1:2) = [protind(jj),...
                            tempresult_dist(ll).ProteinID(2) + tempnewprotbatchindoffset];
                    end
                    tempnewprotbatchindoffset = tempnewprotbatchindoffset + length(tempnewprotbatch_dist);
                    result = [result,tempresult_dist];
                    tempnewprotbatch = [tempnewprotbatch;tempnewprotbatch_dist];
                end
            end
        else
            for kk = 1:length(thisprotbatch)
                thisprotbatch_pepfrag{kk} = fragpep(thisprotbatch{kk},1,'bcyz','');
            end
            if input.dohcdtrigger
                SASSOind = 1:size(SASSO,1);
                [result,tempnewprotbatch] = spmd_searchogly(thisprotbatch,thisprotbatch_pepfrag,...
                    FASTAbatch{jj},thisprotbatchmass,PTM,glycancombi,glycancombimass,scoreoptions,SASSO,colnames,...
                    agpoptions,scannum,spectra,charge,precmz,retime,exptmass,precmass,quant,SASSOind);
            else
                SASSO = ms1searchscannum(:);
                SASSOind = 1:length(ms1searchscannum);
                colnames = {'HCD'};
                [result,tempnewprotbatch] = spmd_searchogly(thisprotbatch,thisprotbatch_pepfrag,...
                    FASTAbatch{jj},thisprotbatchmass,PTM,glycancombi,glycancombimass,scoreoptions,SASSO,colnames,...
                    agpoptions,scannum,spectra,charge,precmz,retime,exptmass,precmass,quant,SASSOind);
            end
        end
        allpepdata.glypep{protind(jj)} = tempnewprotbatch;
    end
    savescoreresults(result,scoreoptions,input,PTM,allpepdata,triggerdata,...
        outputdir,outputfilename,protindse,'save');
end
end