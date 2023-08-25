function isomericglypepresult = calcsgp_legacy(result,sgps,scans,displaydataind,allspectra,allscan,allchg,...
    denoisingoptHCD,denoisingoptCID,denoisingoptEThcD,...
    scoreintdata,allfragments,searchrowind,options)
wtbar = waitbar(0,'Calculating');
sectionsize = 9000;
isomericresult = [];
tempresultsto = cell(sectionsize,29);
% colnames = scoreintdata.colnames;
currentind = 0;  % initialize write index
nummindiagionfragmode = 1;

fldnms = fieldnames(result);
scoreoptions = scoreintdata.scoreoptions;
glylib = scoreintdata.ptmisomers;
new_ptmfragsto = allfragments.ptmfragsto(scoreintdata.scoreoptions.fragnum(:,2),:);
new_ptmseq = allfragments.ptmseq;
scoreintdata.ptmfragsto = allfragments.ptmfragsto(scoreintdata.scoreoptions.fragnum(:,2),:);
scoreintdata.ptmfragstublen = allfragments.ptmfragstublen(scoreintdata.scoreoptions.fragnum(:,2),:);
scoreintdata.ptmmass = allfragments.ptmmass;
scoreintdata.ptmtype = allfragments.ptmtype;
scoreintdata.ptmseq = allfragments.ptmseq;
analyzefragmode = scoreoptions.analyzefragmode;
ms2tol = scoreoptions.ms2tol;
ms2tolunit = scoreoptions.ms2tolunit;

displaydataind = displaydataind(searchrowind,:);
for ii = 1:size(displaydataind,1)
    if mod(ii,100) == 0
        waitbar(ii/size(displaydataind,1),wtbar);
    end
    tempsgp = sgps{displaydataind(ii,1)};
    tempscan = scans(displaydataind(ii,:));
    tempchg = allchg(allscan == tempscan);
    tempspectrum = allspectra{allscan == tempscan};
    tempfragmode = result(displaydataind(ii,1)).Fragmode;
    tempmethodind = ismember(upper(analyzefragmode),upper(tempfragmode));
    tempms2tol = ms2tol(tempmethodind);
    tempms2tolunit = ms2tolunit{tempmethodind};
    [p,g_ori,m] = breakGlyPep(tempsgp);
    glycan = g_ori.struct;
    rowind = find(ismember(glylib(:,1),glycan));
    if ~any(rowind)
        error('Isomer table misalignment');
    end
    isomers = glylib(rowind,:);
    isomers = isomers(~cellfun(@isempty,isomers));
    bestisomers = '';
    outputtoresult = true;
    if length(isomers) > 1
        %% SPECIAL STRUCT FEATURES
        numspecialfeatures = 9;
        structdetail_specialstruct = false(numspecialfeatures,length(isomers));
        sgpseq_specialstruct = cell(numspecialfeatures,1);
        sgpseq_specialstruct{3} = {'{n{f}{h}}','{n{h}{f}}'};
        sgpseq_specialstruct{4} = {'{n{f}{h{s}}}','{n{h{s}}{f}}'};
        sgpseq_specialstruct{7} = {'{n{h{f}{n}}}','{n{h{n}{f}}}'};
        sgpseq_specialstruct{8} = {'{n{h{f}{h}}}','{n{h{h}{f}}}'};
        sgpseq_specialstruct{9} = {'{h{f}}'};
        
        %         structdetail_corefuc
        %         structdetail_bisect
        %         structdetail_LeX
        %         structdetail_sLeX
        %         structdetail_termsia
        %         structdetail_diLacNAc
        %         structdetail_bldgpA
        %         structdetail_bldgpB
        %         structdetail_bldgpO
        for jj = 1:length(isomers)
            bondmap = getglycanbondmap(isomers{jj});
            distance = getdistance(isomers{jj});
            monosac = regexp(isomers{jj},'[hnsf]','match');
            fdist = distance(strcmpi(monosac,'f'));
            ndist = distance(strcmpi(monosac,'n'));
            sdist = distance(strcmpi(monosac,'s'));
            npos = find(strcmpi(monosac,'n'));
            for kk = 1:numspecialfeatures
                switch kk
                    case 1
                        if any(fdist == 2)
                            structdetail_specialstruct(kk,jj) = true;
                        end
                    case 2
                        if any(ndist == 4)
                            structdetail_specialstruct(kk,jj) = true;
                        end
                    case 3
                        for ll = 1:length(sgpseq_specialstruct{kk})
                            if any(strfind(isomers{jj},sgpseq_specialstruct{kk}{ll}))
                                structdetail_specialstruct(kk,jj) = true;
                            end
                        end
                    case 4
                        for ll = 1:length(sgpseq_specialstruct{kk})
                            if any(strfind(isomers{jj},sgpseq_specialstruct{kk}{ll}))
                                structdetail_specialstruct(kk,jj) = true;
                            end
                        end
                    case 5
                        if any(sdist > 1)
                            structdetail_specialstruct(kk,jj) = true;
                        end
                    case 6
                        tempnpos = setdiff(npos,1);
                        for ll = 1:length(tempnpos)
                            if any(ismember('n',monosac(logical(bondmap(tempnpos(ll),:)))))
                                structdetail_specialstruct(kk,jj) = true;
                            end
                        end
                    case 7
                        for ll = 1:length(sgpseq_specialstruct{kk})
                            if any(strfind(isomers{jj},sgpseq_specialstruct{kk}{ll}))
                                structdetail_specialstruct(kk,jj) = true;
                            end
                        end
                    case 8
                        for ll = 1:length(sgpseq_specialstruct{kk})
                            if any(strfind(isomers{jj},sgpseq_specialstruct{kk}{ll}))
                                structdetail_specialstruct(kk,jj) = true;
                            end
                        end
                    case 9
                        for ll = 1:length(sgpseq_specialstruct{kk})
                            if any(strfind(isomers{jj},sgpseq_specialstruct{kk}{ll}))
                                structdetail_specialstruct(kk,jj) = true;
                            end
                        end
                end
            end
        end
        
        %% PREPARE SPECTRUM
        if strcmpi(tempfragmode,'HCD')
            tempspectrum = spectradenoising(tempspectrum,{'HCD'},tempchg,...
                glypepMW(tempsgp),...
                tempms2tol,tempms2tolunit,denoisingoptHCD);
            theofrag_HCDsto = cell(length(isomers),1);
        elseif strcmpi(tempfragmode,'CID')
            tempspectrum = spectradenoising(tempspectrum,{'CID'},tempchg,...
                glypepMW(tempsgp),tempms2tol,tempms2tolunit,denoisingoptCID);
            theofrag_CIDsto = cell(length(isomers),1);
        elseif strcmpi(tempfragmode,'ETHCD')
            tempspectrum = spectradenoising(tempspectrum,{'EThcD'},tempchg,...
                glypepMW(tempsgp),tempms2tol,tempms2tolunit,denoisingoptEThcD);
            theofrag_EThcD_cofragsto = cell(length(isomers),1);
            theofrag_EThcD_glyonlysto = cell(length(isomers),1);
        end
        
        %% PREPARE THEO FRAG ION
        isomerProteinIDsto = cell(length(isomers),1);
        isomerfrags = [];
        for jj = 1:length(isomers)
            thisisomerind = find(ismember(new_ptmseq,isomers{jj}));
            newprotid = str2num(result(displaydataind(ii,1)).ProteinID);
            newprotid1 = newprotid(1:2);
            newprotid2 = newprotid(3:end);
            for kk = 1:length(newprotid2)
                currentptmid = newprotid2(kk);
                if currentptmid == rowind  % this is the one to swap
                    newprotid2(kk) = thisisomerind;
                else
                    newprotid2(kk) = find(ismember(new_ptmseq,glylib{newprotid2(kk)}));
                end
            end
            newprotid = [newprotid1,newprotid2];
            isomerProteinIDsto{jj} = newprotid;
            if isempty(new_ptmfragsto{1,thisisomerind})
                [~,tempptmmass,~,...
                    tempptmtype,tempptmfragsto,tempptmfragstublen] = ...
                    theoptmfrag(isomers(jj),scoreoptions.fragnum(:,2),scoreoptions.analyzefragmode,...
                    scoreoptions);
                scoreintdata.ptmfragsto(:,thisisomerind) = tempptmfragsto;
                scoreintdata.ptmfragstublen(:,thisisomerind) = tempptmfragstublen;
                scoreintdata.ptmmass(thisisomerind) = tempptmmass;
                scoreintdata.ptmtype(thisisomerind) = tempptmtype;
            end
            temptheofrag = createtheofragsto({newprotid},{tempfragmode},scoreintdata);
            switch upper(tempfragmode)
                case 'HCD'
                    theofrag_HCD = temptheofrag{1};
                    for kk = 1:length(theofrag_HCD)
                        theofrag_HCD(kk).isomerind = jj;
                    end
                    theofrag_HCDsto{jj} = theofrag_HCD;
                case 'CID'
                    theofrag_CID = temptheofrag{1};
                    for kk = 1:length(theofrag_CID)
                        theofrag_CID(kk).isomerind = jj;
                    end
                    theofrag_CIDsto{jj} = theofrag_CID;
                case 'ETHCD'
                    theofrag_EThcD_cofrag = temptheofrag{1};
                    theofrag_EThcD_glyonly = theofrag_EThcD_cofrag;
                    theofrag_EThcD_glyonly = theofrag_EThcD_glyonly([theofrag_EThcD_glyonly.ngFrag] > 0 & ...
                        [theofrag_EThcD_glyonly.npFrag] == 0);
                    for kk = 1:length(theofrag_EThcD_cofrag)
                        theofrag_EThcD_cofrag(kk).isomerind = jj;
                    end
                    for kk = 1:length(theofrag_EThcD_glyonly)
                        theofrag_EThcD_glyonly(kk).isomerind = jj;
                    end
                    theofrag_EThcD_cofragsto{jj} = theofrag_EThcD_cofrag;
                    theofrag_EThcD_glyonlysto{jj} = theofrag_EThcD_glyonly;
            end
        end
        
        %% DELETE BAD ISOMERS - DIAGNOSTIC ION BASED
        diagnosticmass_B = [];
        diagnosticmasstyp_B = [];
        diagnosticmass_Y = [];
        diagnosticmasstyp_Y = [];
        diagnosticmass_Blacdinac = [];
        diagnosticmasstyp_Blacdinac = [];
        diagnosticmass_Ylacdinac = [];
        diagnosticmasstyp_Ylacdinac = [];
        if any(structdetail_specialstruct(1,:))
            g_diag = g_ori;
            g_diag.struct = '{n{f}}';
            diag_1 = glypepMW(joinGlyPep(p,g_diag,m)) + 1.007825032;
            g_diag.struct = '{n{f}{n}}';
            diag_2 = glypepMW(joinGlyPep(p,g_diag,m)) + 1.007825032;
            diagnosticmass_Y = [diagnosticmass_Y;diag_1;diag_2];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;1;1];
        end
        if any(structdetail_specialstruct(2,:))
            g_diag = g_ori;
            g_diag.struct = '{n{n{h{n}}}}';
            diag_1 = glypepMW(joinGlyPep(p,g_diag,m)) + 1.007825032;
            g_diag.struct = '{n{f}{n{h{n}}}}';
            diag_2 = glypepMW(joinGlyPep(p,g_diag,m)) + 1.007825032;
            diagnosticmass_Y = [diagnosticmass_Y;diag_1;diag_2];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;2;2];
        end
        if any(structdetail_specialstruct(3,:))
            diag_1 = glyMW('{n{f}{h}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;3];
        end
        if any(structdetail_specialstruct(4,:))
            diag_1 = glyMW('{n{f}{h{s}}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;4];
        end
        if any(structdetail_specialstruct(5,:))
            diag_1 = glyMW('{s}') + 1.007825032 - 18.0105647;
            diag_2 = glyMW('{h{s}}') + 1.007825032 - 18.0105647;
            diag_3 = glyMW('{n{h{s}}}') + 1.007825032 - 18.0105647;
            diag_4 = glyMW('{n{n{s}}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1;diag_2;diag_3;diag_4];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;5;5;5;5];
        end
        if any(structdetail_specialstruct(6,:))
            tempdiagnosticmass_Blacdinac = [];
            tempdiagnosticmass_Ylacdinac = [];
            for jj = 1:length(isomers)
                bondmap = getglycanbondmap(isomers{jj});
                monosac = regexp(isomers{jj},'[hnsf]','match');
                ndist = distance(strcmpi(monosac,'n'));
                antennarootNind = find(ndist == 5);
                for kk = 1:length(antennarootNind)
                    Nchildren = glytreetracker(bondmap,antennarootNind(kk),[],'down');
                    if ismember('n',monosac(Nchildren))
                        tempdiagnosticmass_Blacdinac = [tempdiagnosticmass_Blacdinac;...
                            glyMW(strjoin(monosac(Nchildren))) + 1.007825032 - 18.0105647];
                        tempdiagnosticmass_Ylacdinac = [tempdiagnosticmass_Ylacdinac;...
                            glyMW(isomers{jj}) - glyMW(strjoin(monosac(Nchildren))) + 1.007825032 + 18.0105647];
                    end
                end
            end
            diagnosticmass_Blacdinac = [diagnosticmass_Blacdinac;tempdiagnosticmass_Blacdinac];
            diagnosticmasstyp_Blacdinac = [diagnosticmasstyp_Blacdinac;ones(size(tempdiagnosticmass_Blacdinac)) * 6];
            diagnosticmass_Ylacdinac = [diagnosticmass_Ylacdinac;tempdiagnosticmass_Ylacdinac];
            diagnosticmasstyp_Ylacdinac = [diagnosticmasstyp_Ylacdinac;ones(size(tempdiagnosticmass_Ylacdinac)) * 6];
        end
        if any(structdetail_specialstruct(7,:))
            diag_1 = glyMW('{n{h{f}{n}}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;7];
        end
        if any(structdetail_specialstruct(8,:))
            diag_1 = glyMW('{h{f}{h}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;8];
        end
        if any(structdetail_specialstruct(9,:))
            diag_1 = glyMW('{h{f}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;9];
        end
        if ~isempty(diagnosticmass_B)
            found_specialstruct = false(numspecialfeatures,1);
            diagmatch_B = calcithscore(tempspectrum,diagnosticmass_B,1,...
                tempms2tol,tempms2tolunit,3,struct('maxlag',1,'selectpeak',[]));
            diagmatch_Y = calcithscore(tempspectrum,diagnosticmass_Y,tempchg,...
                tempms2tol,tempms2tolunit,3,struct('maxlag',1,'selectpeak',[]));
            diagmatch_Blacdinac = calcithscore(tempspectrum,diagnosticmass_Blacdinac,1,...
                tempms2tol,tempms2tolunit,3,struct('maxlag',1,'selectpeak',[]));
            diagmatch_Ylacdinac = calcithscore(tempspectrum,diagnosticmass_Ylacdinac,1,...
                tempms2tol,tempms2tolunit,3,struct('maxlag',1,'selectpeak',[]));
            matcheddiagiontyp = [diagnosticmasstyp_B(logical(diagmatch_B.ionmatchindex));...
                diagnosticmasstyp_Y(logical(diagmatch_Y.ionmatchindex))];
            matcheddiagiontyp_Blacdinac = logical(diagmatch_Blacdinac.ionmatchindex);
            matcheddiagiontyp_Ylacdinac = logical(diagmatch_Ylacdinac.ionmatchindex);
            for jj = 1:numspecialfeatures
                if jj == 6
                    featurefoundbydiagion = matcheddiagiontyp_Blacdinac & matcheddiagiontyp_Ylacdinac;
                    if any(featurefoundbydiagion)  % found feature 1~9
                        found_specialstruct(jj) = true;
                    end
                elseif any(structdetail_specialstruct(jj,:) == 1)
                    % mixed - isomer may or maynot has struct feature
                    featurefoundbydiagion = matcheddiagiontyp == jj;
                    if any(featurefoundbydiagion)  % found feature 1~9
                        found_specialstruct(jj) = true;
                    end
                end
            end
        else
            found_specialstruct = true(numspecialfeatures,1);
        end
        isomerkeepind = false(size(isomers));
        isomerhasstructfeatures = any(any(structdetail_specialstruct));
        if isomerhasstructfeatures  % any isomer has structfeatures
            for jj = 1:length(isomers)
                if any(structdetail_specialstruct(:,jj))  % this isomer has some structfeatures
                    keepthisisomer = ~xor(structdetail_specialstruct(:,jj),found_specialstruct >= nummindiagionfragmode);  % feature found
                    if ~any(~keepthisisomer)
                        isomerkeepind(jj) = true;
                    end
                end
            end
            if ~any(isomerkeepind)
                isomerkeepind = ~xor(any(structdetail_specialstruct),isomerkeepind);  % 1st half: each isomer has sf or not, 2nd half: all false
                % above operation keeps only those without any
                % structfeatures
            end
        else  % no isomer has structfeatures
            isomerkeepind = true(size(isomers));
        end
        isomers = isomers(isomerkeepind);
        
        %% ISOMER RANKING
        if isempty(isomers)
            outputtoresult = false;
        else
            switch upper(tempfragmode)
                case 'HCD'
                    theofrag_HCDsto = theofrag_HCDsto(isomerkeepind);
                    for jj = 1:length(isomers)
                        [theofrag_HCDsto{jj}.isomerind] = deal(jj);
                        isomerfrags = [isomerfrags,theofrag_HCDsto{jj}];
                    end
                case 'CID'
                    theofrag_CIDsto = theofrag_CIDsto(isomerkeepind);
                    for jj = 1:length(isomers)
                        [theofrag_CIDsto{jj}.isomerind] = deal(jj);
                        isomerfrags = [isomerfrags,theofrag_CIDsto{jj}];
                    end
                case 'ETHCD'
                    theofrag_EThcD_cofragsto = theofrag_EThcD_cofragsto(isomerkeepind);
                    theofrag_EThcD_glyonlysto = theofrag_EThcD_glyonlysto(isomerkeepind);
                    for jj = 1:length(isomers)
                        [theofrag_EThcD_cofragsto{jj}.isomerind] = deal(jj);
                        [theofrag_EThcD_glyonlysto{jj}.isomerind] = deal(jj);
                        isomerfrags = [isomerfrags,theofrag_EThcD_glyonlysto{jj}];
                    end
            end
            isomerProteinIDsto = isomerProteinIDsto(isomerkeepind);
            
            %% MATCH EXPT SPECTRUM
            isomerfrags_mz = [];
            isomerfrags_mzisomerind = [];
            isomerfrags_chg = [];
            isomerfrags_fragind = [];
            isomerfrags_originalisomerind = [isomerfrags.isomerind];
            for jj = 1:tempchg  % get the biglist
                tempisomerfrags_mz = ([isomerfrags.mz] - 1.007825032)/jj + 1.007825032;
                isomerfrags_mz = [isomerfrags_mz;tempisomerfrags_mz(:)];
                isomerfrags_mzisomerind = [isomerfrags_mzisomerind;[isomerfrags.isomerind]];
                isomerfrags_chg = [isomerfrags_chg;ones(length(isomerfrags),1) * jj];
                isomerfrags_fragind = [isomerfrags_fragind;(1:length(isomerfrags))'];
            end
            [~,~,isomerfrags_uniqueness] = unique(isomerfrags_mz,'stable');
            match = calcithscore(tempspectrum,[isomerfrags.mz],tempchg,...
                tempms2tol,tempms2tolunit,2,...
                struct('maxlag',scoreoptions.maxlag,'selectpeak',[]));
            medianpeak = median(tempspectrum(:,2));
            matchedpeakind = match.peakmatchindex;
            matchedpeakind_type = zeros(size(matchedpeakind,1),1);  % temp variable, store info for potential use
            uniisomers_numerical = zeros(1,length(isomers));
            sharedisomers_numerical = zeros(1,length(isomers));
            uniisomers_quant = zeros(1,length(isomers));
            sharedisomers_quant = zeros(1,length(isomers));
            allmatchedpeakind = [];
            allmatchedpeakind_whichpeak = [];
            for jj = 1:size(matchedpeakind,1)
                if ~isempty(matchedpeakind{jj})
                    allmatchedpeakind = [allmatchedpeakind;matchedpeakind{jj}];
                    allmatchedpeakind_whichpeak = [allmatchedpeakind_whichpeak;...
                        repmat(jj,size(matchedpeakind{jj},1),1)];
                end
            end
            [allmatchedpeakind,~,allmatchedpeakind_uniind] = unique(allmatchedpeakind,'rows');
            tempimportantstructfeaturesfound_isomerid = cell(size(allmatchedpeakind,1),1);
            for jj = 1:size(allmatchedpeakind,1)  % fragind, chg
                tempmatchedunifragind = allmatchedpeakind(jj,1);  % matched frag ind
                tempmatchedunifragchg = allmatchedpeakind(jj,2);  % maeched frag chg
                tempmatchedpeakind = unique(allmatchedpeakind_whichpeak(allmatchedpeakind_uniind == jj));
                tempmatchedunifragind_inbig = isomerfrags_fragind == tempmatchedunifragind & ...
                    isomerfrags_chg == tempmatchedunifragchg;  % find its absolute position in biglist
                tempmatchedfraguniqueind = isomerfrags_uniqueness(tempmatchedunifragind_inbig);
                matchedfragind_fororiginal = isomerfrags_fragind(ismember(isomerfrags_uniqueness,tempmatchedfraguniqueind));
                matchedorifrag  = isomerfrags(matchedfragind_fororiginal);
                matchedorifragchg = isomerfrags_chg(ismember(isomerfrags_uniqueness,tempmatchedfraguniqueind));
                keepthisfrag = true(size(matchedorifrag));
                sigionquant = sum(log2(tempspectrum(tempmatchedpeakind,2)/medianpeak));
                for kk = 1:length(matchedorifrag)
                    temptype = matchedorifrag(kk).type;
                    if strcmpi(temptype(1:2),'-b')  % Bion
                        if matchedorifragchg(kk) > 1
                            keepthisfrag(kk) = false;
                        end
                    elseif strcmpi(temptype(1:2),'-y')  % Bion
                        if matchedorifragchg(kk) > 3
                            keepthisfrag(kk) = false;
                        end
                    end
                end
                matchedfragind_fororiginal = matchedfragind_fororiginal(keepthisfrag);
                [tempmatchedisomerind,~,tempmatchedisomerindind2] = ...
                    unique(isomerfrags_originalisomerind(matchedfragind_fororiginal));
                if length(tempmatchedisomerind) == 1  % unique
                    matchedpeakind_type(jj) = 1;
                    uniisomers_numerical(tempmatchedisomerind) = uniisomers_numerical(tempmatchedisomerind) + ...
                        length(tempmatchedisomerindind2);
                    uniisomers_quant(tempmatchedisomerind) = uniisomers_quant(tempmatchedisomerind) + ...
                        sigionquant;
                    tempimportantstructfeaturesfound_isomerid{jj} = tempmatchedisomerind;
                elseif length(tempmatchedisomerind) == length(isomers)  % common
                    matchedpeakind_type(jj) = -1;
                elseif ~isempty(tempmatchedisomerind)  % shared
                    matchedpeakind_type(jj) = 2;
                    tempisomerind = [];
                    for kk = 1:length(tempmatchedisomerind)
                        numisomerfound = sum(tempmatchedisomerindind2 == tempmatchedisomerindind2(kk));
                        sharedisomers_numerical(tempmatchedisomerind(kk)) = ...
                            sharedisomers_numerical(tempmatchedisomerind(kk)) + numisomerfound;
                        sharedisomers_quant(tempmatchedisomerind(kk)) = sharedisomers_quant(tempmatchedisomerind(kk)) + ...
                            sigionquant;
                        tempisomerind = [tempisomerind,tempmatchedisomerind(kk)];
                    end
                end
            end
            % Merge - If EThcD_cofrag is needed
            tempbestisomers_string = '';
            tempisomers_noncommonfound_numerical = uniisomers_numerical + sharedisomers_numerical;
            tempisomers_numnoncommonfragfound_quant = uniisomers_quant + sharedisomers_quant;
            if any(tempisomers_noncommonfound_numerical)
                foundisomersind = (tempisomers_noncommonfound_numerical > 0);
                foundisomers = isomers(foundisomersind);
                foundisomers_numnoncommonfragfound_quant = tempisomers_numnoncommonfragfound_quant(foundisomersind);
                [foundisomers_numnoncommonfragfound_quant,ind] = sort(foundisomers_numnoncommonfragfound_quant,'descend');
                foundisomers = foundisomers(ind);
                if length(foundisomers_numnoncommonfragfound_quant) == 1
                    tempbestisomers = foundisomers;
                elseif any(foundisomers_numnoncommonfragfound_quant < 0.8 * max(foundisomers_numnoncommonfragfound_quant))
                    [idx,~] = kmeans(foundisomers_numnoncommonfragfound_quant(:),2);
                    tempbestisomers = foundisomers(idx == 2);
                else
                    tempbestisomers = foundisomers;
                end
            else
                tempbestisomers = isomers;
            end
            for jj = 1:length(tempbestisomers)
                tempbestisomers_string = [tempbestisomers_string,tempbestisomers{jj},', '];
            end
            bestisomers = tempbestisomers_string(1:end - 2);
        end
    end
    %% Isomeric glypep scoring
    if outputtoresult
        if isempty(bestisomers)
            currentind = currentind + 1;
            thisisomerind = find(ismember(new_ptmseq,isomers{1}));
            newprotid = str2num(result(displaydataind(ii,1)).ProteinID);
            newprotid1 = newprotid(1:2);
            newprotid2 = newprotid(3:end);
            for kk = 1:length(newprotid2)
                currentptmid = newprotid2(kk);
                if currentptmid == rowind  % this is the one to swap
                    newprotid2(kk) = thisisomerind;
                end
            end
            newprotid = [newprotid1,newprotid2];
            tempresult = result(displaydataind(ii,1));
            if strcmpi(tempresult.Fragmode,'ethcd')
                tempresult.Fragmode = 'EThcD';
            end
            tempresult.ProteinID = newprotid;
            if currentind > sectionsize
                isomericresult = [isomericresult;tempresultsto];
                tempresultsto = cell(sectionsize,29);
                currentind = 1;
            end
            for kk = 1:length(fldnms)
                tempresultsto{currentind,kk} = tempresult.(fldnms{kk});
            end
        else
            switch upper(tempfragmode)
                case 'HCD'
                    theofragstoptmseq = cell(1,length(theofrag_HCDsto));
                    for jj = 1:length(theofragstoptmseq)
                        theofragstoptmseq{jj} = theofrag_HCDsto{jj}(1).original;
                    end
                case 'CID'
                    theofragstoptmseq = cell(1,length(theofrag_CIDsto));
                    for jj = 1:length(theofrag_CIDsto)
                        theofragstoptmseq{jj} = theofrag_CIDsto{jj}(1).original;
                    end
                case 'ETHCD'
                    theofragstoptmseq = cell(1,length(theofrag_EThcD_cofragsto));
                    for jj = 1:length(theofrag_EThcD_cofragsto)
                        theofragstoptmseq{jj} = theofrag_EThcD_cofragsto{jj}(1).original;
                    end
            end
            tempbestisomers = strsplit(bestisomers,', ');
            newisomericglycan = setdiff(tempbestisomers,g_ori.struct);
            tempresult = result(displaydataind(ii,1));
            tempresult.ProteinID = isomerProteinIDsto{1};
            if ismember(g_ori.struct,tempbestisomers)
                currentind = currentind + 1;
                if currentind > sectionsize
                    isomericresult = [isomericresult;tempresultsto];
                    tempresultsto = cell(sectionsize,29);
                    currentind = 1;
                end
                for kk = 1:length(fldnms)
                    tempresultsto{currentind,kk} = tempresult.(fldnms{kk});
                end
            end
            if ~isempty(newisomericglycan)
                for jj = 1:length(newisomericglycan)
                    thisnewisomericglycan = newisomericglycan{jj};
                    g_alt = g_ori;
                    g_alt.struct = thisnewisomericglycan;
                    thisnewglypepseq = joinGlyPep(p,g_alt,m);
                    sgpind = ismember(theofragstoptmseq,thisnewglypepseq);
                    isomericglypepprotid = isomerProteinIDsto{sgpind};
                    switch upper(tempfragmode)
                        case 'HCD'
                            theofrag = theofrag_HCDsto{sgpind};
                        case 'CID'
                            theofrag = theofrag_CIDsto{sgpind};
                        case 'ETHCD'
                            theofrag = theofrag_EThcD_cofragsto{sgpind};
                    end
                    [spectrascores,~] = match1by1(thisnewglypepseq,tempresult.Theo,...
                        theofrag,isomericglypepprotid,tempscan,...
                        tempresult.Expt,tempresult.Mono,tempresult.Retime,tempspectrum,tempchg,...
                        tempms2tol,tempms2tolunit,tempfragmode,scoreintdata.sliminput.fragnum,tempresult.MDiff,...
                        tempresult.Quant,tempresult.Protein,struct('maxlag',scoreoptions.maxlag,'selectpeak',[],'proceed_override',false),false);
                    spectrascores = rmfield(spectrascores,'Pepcov');
                    %added fldnms_new
                    fldnms_new = fieldnames(result);
                    currentind = currentind + 1;
                    if currentind > sectionsize
                        isomericresult = [isomericresult;tempresultsto];
                        tempresultsto = cell(sectionsize,29);
                        currentind = 1;
                    end
                    for kk = 1:length(fldnms_new)
                        tempresultsto{currentind,kk} = spectrascores.(fldnms_new{kk});
                    end
                end
            end
        end
    else
%         tempresultsto = {};
    end
end
isomericresult = [isomericresult;tempresultsto];
if ~isempty(isomericresult)
    isomericresult = isomericresult(~cellfun(@isempty,isomericresult(:,1)),:);
    ProteinIDcol = ismember(fldnms,'ProteinID');
    for ii = 1:size(isomericresult,1)
        tempProteinID = isomericresult{ii,ProteinIDcol};
        if isnumeric(tempProteinID)
            isomericresult{ii,ProteinIDcol} = num2str(tempProteinID);
        end
    end
    isomericglypepresult = cell2struct(isomericresult,fldnms,2)';
else
    isomericglypepresult = [];
end
delete(wtbar);
end