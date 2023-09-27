function isomericglypepresult = calcsgp(result,sgps,scans,displaydataind,allspectra,allscan,allchg,...
    denoisingoptHCD,denoisingoptCID,denoisingoptETD,denoisingoptETciD,denoisingoptEThcD,...
    scoreintdata,allfragments,searchrowind,options)
% 'result' has filtered output from 1_1 file that contains only
% combiES>0.4, and addition filter settings in GlycoPATConstant.m
% 'scans' contains the HCD scans that are paired with CID and ETHCD using
% SASSO
% 'sgp' has SGP1.0 corresponding to the individual scans
% wtbar = waitbar(0,'Calculating');
ct=0;
sectionsize = 9000;
isomericresult = [];
tempresult = cell(sectionsize,30);
colnames = scoreintdata.colnames;
currentind = 1 - length(colnames);  % initialize write index
nummindiagionfragmode = min(GlycoPATConstant.fbi_nummmindiagionfragmode,...
    length(scoreintdata.scoreoptions.analyzefragmode));
diagionlogicandor = GlycoPATConstant.fbi_diagnionlogic;  % 1: AND; 2: OR

fldnms = fieldnames(result);
scoreoptions = scoreintdata.scoreoptions;
glylib = scoreintdata.ptmisomers;
new_ptmfragsto = allfragments.ptmfragsto;
new_ptmseq = allfragments.ptmseq;
scoreintdata.ptmfragsto = allfragments.ptmfragsto;
scoreintdata.ptmfragstublen = allfragments.ptmfragstublen;
scoreintdata.ptmmass = allfragments.ptmmass;
scoreintdata.ptmtype = allfragments.ptmtype;
scoreintdata.ptmseq = allfragments.ptmseq;
analyzefragmode = scoreoptions.analyzefragmode;
ms2tol = scoreoptions.ms2tol;
ms2tolunit = scoreoptions.ms2tolunit;
methodind_CID = ismember(upper(analyzefragmode),'CID');
methodind_HCD = ismember(upper(analyzefragmode),'HCD');
methodind_ETD = ismember(upper(analyzefragmode),'ETD');
methodind_ETciD = ismember(upper(analyzefragmode),'ETCID');
methodind_EThcD = ismember(upper(analyzefragmode),'ETHCD');

displaydataind = displaydataind(searchrowind,:);            % scan numbers triggered by single precursor corresponding to each candidate hit, typicaly a SASSO row
displaydataind_colCID = ismember(upper(colnames),'CID');
displaydataind_colHCD = ismember(upper(colnames),'HCD');
displaydataind_colETD = ismember(upper(colnames),'ETD');
displaydataind_colETciD = ismember(upper(colnames),'ETCID');
displaydataind_colEThcD = ismember(upper(colnames),'ETHCD');
for ii = 1:size(displaydataind,1) % 172 groups of scans
    %     if mod(ii,100) == 0
    %         waitbar(ii/size(displaydataind,1),wtbar);
    %     end
    tempsgp = sgps{displaydataind(ii,1)};                   % sgp identified in the first set of runs with unique glycans, i.e. 1_1 run
    tempscans = scans(displaydataind(ii,:));
    gpmw = result(displaydataind(ii,1)).Theo;
    [p,g_ori,m] = breakGlyPep(tempsgp);
    glycan = g_ori.struct;                                  % candidate_glycan identied in 1_1 corresponding to displaydataind
    rowind = find(ismember(glylib(:,1),glycan));
    if ~any(rowind)
        error('Isomer table misalignment');
    end
    isomers = glylib(rowind,:);
    isomers = isomers(~cellfun(@isempty,isomers));          % identifies isomers corresponding to candidate_glycan
    bestisomers = '';
    outputtoresult = true;
    if length(isomers) > 1                                  % which of the isomers contain the special structural features that are present in the experimental spectrum
        %% SPECIAL STRUCT FEATURES
        numspecialfeatures = 13;
        structdetail_specialstruct = false(numspecialfeatures,length(isomers));
        sgpseq_specialstruct = cell(numspecialfeatures,1);
        %         structdetail_corefuc (sgpseq_specialstruct{1})
        %         structdetail_bisect  (sgpseq_specialstruct{2})
        sgpseq_specialstruct{3} = {'{n{f}{h}}','{n{h}{f}}'};  % Lewis-X
        %sgpseq_specialstruct{4} = {'{n{f}{h{s}}}','{n{h{s}}{f}}','{n{f}}','{n{h{s}}}'}; % sialyl Lewis-X
        sgpseq_specialstruct{4} = {'{n{f}{h{s}}}','{n{h{s}}{f}}'}; % sialyl Lewis-X
        sgpseq_specialstruct{5} = {'{n{h{s}}}','{s}'};        % terminal sialic acid
        sgpseq_specialstruct{6} = {'{n{n}}'};                 % diLacNAc
        sgpseq_specialstruct{7} = {'{h{f}{n}}','{h{n}{f}}'};  % bldgpA
        sgpseq_specialstruct{8} = {'{h{f}{h}}','{h{h}{f}}'};  % bldgpB
        sgpseq_specialstruct{9} = {'{h{f}}'};           % blood gp.
        sgpseq_specialstruct{10} = {'{n{f}{h{f}}}'};    % Lewis Y/b
        sgpseq_specialstruct{11} = {'{n{f}{n}}','{n{n}{f}}'};   % fucosylated diLacNAc
        sgpseq_specialstruct{12} = {'{n{n{s}}}'};               % sialylated diLacNAc
        sgpseq_specialstruct{13} = {'{h{h}}','{h{h{h}}}'};      % hybrid

        
        for jj = 1:length(isomers)                      % This is to check if the candidate isomers sgp structure has diagnostic ions. If any exist then it is marked true, else false
            bondmap = getglycanbondmap(isomers{jj});    % get connection table for glycan structure
            distance = getdistance(isomers{jj});        % returns the distance from monosaccharide to reducing end
            monosac = regexp(isomers{jj},'[hnsf]','match');
            fdist = distance(strcmpi(monosac,'f'));
            ndist = distance(strcmpi(monosac,'n'));
            sdist = distance(strcmpi(monosac,'s'));
            npos = find(strcmpi(monosac,'n'));
            for kk = 1:numspecialfeatures
                switch kk
                    case 1
                        if any(fdist == 2)    % core fucose
                            structdetail_specialstruct(kk,jj) = true;
                        end
                    case 2
                        if any(ndist == 4)    % bisecting
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
                        for ll = 1:length(sgpseq_specialstruct{kk})
                            if any(strfind(isomers{jj},sgpseq_specialstruct{kk}{ll}))
                                structdetail_specialstruct(kk,jj) = true;
                            end
                        end
                    case 6
                        for ll = 1:length(sgpseq_specialstruct{kk})
                            if any(strfind(isomers{jj},sgpseq_specialstruct{kk}{ll}))
                                structdetail_specialstruct(kk,jj) = true;
                            end
                        end
                        %                         tempnpos = setdiff(npos,1);
                        %                         for ll = 1:length(tempnpos)
                        %                             if any(ismember('n',monosac(logical(bondmap(tempnpos(ll),:)))))
                        %                                 structdetail_specialstruct(kk,jj) = true;
                        %                             end
                        %                         end
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
                    case 10
                        for ll = 1:length(sgpseq_specialstruct{kk})
                            if any(strfind(isomers{jj},sgpseq_specialstruct{kk}{ll}))
                                structdetail_specialstruct(kk,jj) = true;
                            end
                        end
                    case 11
                        for ll = 1:length(sgpseq_specialstruct{kk})
                            if any(strfind(isomers{jj},sgpseq_specialstruct{kk}{ll}))
                                structdetail_specialstruct(kk,jj) = true;
                            end
                        end
                    case 12
                        for ll = 1:length(sgpseq_specialstruct{kk})
                            if any(strfind(isomers{jj},sgpseq_specialstruct{kk}{ll}))
                                structdetail_specialstruct(kk,jj) = true;
                            end
                        end
                    case 13
                        for ll = 1:length(sgpseq_specialstruct{kk})
                            if any(strfind(isomers{jj},sgpseq_specialstruct{kk}{ll}))
                                structdetail_specialstruct(kk,jj) = true;
                            end
                        end
                end
            end
        end

        %% PREPARE SPECTRUM
        if any(methodind_HCD)
            scanHCD = tempscans(displaydataind_colHCD);
            chgHCD = allchg(allscan == scanHCD);
            spectrumHCD = allspectra{allscan == scanHCD};
            spectrumHCD = spectradenoising(spectrumHCD,{'HCD'},chgHCD,...
                glypepMW(tempsgp),...
                ms2tol(methodind_HCD),ms2tolunit{methodind_HCD},denoisingoptHCD);
        end
        if any(methodind_CID)
            scanCID = tempscans(displaydataind_colCID);
            chgCID = allchg(allscan == scanCID);
            spectrumCID = allspectra{allscan == scanCID};
            spectrumCID = spectradenoising(spectrumCID,{'CID'},chgCID,...
                glypepMW(tempsgp),ms2tol(methodind_CID),ms2tolunit{methodind_CID},denoisingoptCID);
        end
        if any(methodind_ETD)
            scanETD = tempscans(displaydataind_colETD);
            chgETD = allchg(allscan == scanETD);
            spectrumETD = allspectra{allscan == scanETD};
            spectrumETD = spectradenoising(spectrumETD,{'ETD'},chgETD,...
                glypepMW(tempsgp),ms2tol(methodind_ETD),ms2tolunit{methodind_ETD},denoisingoptETD);
        end
        if any(methodind_ETciD)
            scanETciD = tempscans(displaydataind_colETciD);
            chgETciD = allchg(allscan == scanETciD);
            spectrumETciD = allspectra{allscan == scanETciD};
            spectrumETciD = spectradenoising(spectrumETciD,{'ETCID'},chgETciD,...
                glypepMW(tempsgp),ms2tol(methodind_ETciD),ms2tolunit{methodind_ETciD},denoisingoptETciD);
        end
        if any(methodind_EThcD)
            scanEThcD = tempscans(displaydataind_colEThcD);
            chgEThcD = allchg(allscan == scanEThcD);
            spectrumEThcD = allspectra{allscan == scanEThcD};
            spectrumEThcD = spectradenoising(spectrumEThcD,{'EThcD'},chgEThcD,...
                glypepMW(tempsgp),ms2tol(methodind_EThcD),ms2tolunit{methodind_EThcD},denoisingoptEThcD);
        end

        %% PREPARE THEO FRAG ION
        theofrag_EThcD_cofragsto = cell(length(isomers),1);         % one theoretical fragmentation list for each isomer and for each fragmentation mode
        theofrag_EThcD_glyonlysto = cell(length(isomers),1);
        theofrag_CIDsto = cell(length(isomers),1);
        theofrag_HCDsto = cell(length(isomers),1);
        theofrag_ETDsto = cell(length(isomers),1);
        theofrag_ETciDsto = cell(length(isomers),1);
        isomerProteinIDsto = cell(length(isomers),1);
        isomerfrags_HCD = [];
        isomerfrags_EThcD_glyonly = [];
        isomerfrags_CID = [];
        isomerfrags_ETD = [];
        isomerfrags_ETciD = [];
        currentglyisomerid = find(ismember(new_ptmseq,isomers{1}));
        for jj = 1:length(isomers)
            thisisomerind = find(ismember(new_ptmseq,isomers{jj}));
            newprotid = str2num(result(displaydataind(ii,1)).ProteinID);
            newprotid1 = newprotid(1:2);        % protein and peptide id
            newprotid2 = newprotid(3:end);      % glycan ID
            for kk = 1:length(newprotid2)
                if newprotid2(kk) == currentglyisomerid
                    newprotid2(kk) = thisisomerind;
                    %                 else
                    %                     newprotid2(kk) = find(ismember(new_ptmseq,scoreintdata.ptmisomers{newprotid2(kk),1}));
                end
            end
            newprotid = [newprotid1,newprotid2];        % protein, peptide id, glycan ID
            isomerProteinIDsto{jj} = newprotid;
            if isempty(new_ptmfragsto{1,thisisomerind})
                [~,tempptmmass,~,...
                    tempptmtype,tempptmfragsto,tempptmfragstublen] = ...
                    theoptmfrag(isomers(jj),scoreoptions.fragnum(:,2),scoreoptions.analyzefragmode,...
                    scoreoptions);                                                      % fragment all PTM appeared in digested protein file
                for kk = 1:size(scoreoptions.fragnum,1)
                    scoreintdata.ptmfragsto(kk,thisisomerind) = tempptmfragsto{kk,1};
                    scoreintdata.ptmfragstublen(kk,thisisomerind) = tempptmfragstublen{kk,1};
                end
                scoreintdata.ptmmass(thisisomerind) = tempptmmass;
                scoreintdata.ptmtype(thisisomerind) = tempptmtype;
            end
            temptheofrag = createtheofragsto({newprotid},colnames,scoreintdata);        % create theoretical fragments for input glycopeptides
            if any(displaydataind_colHCD)
                theofrag_HCD = temptheofrag{1,displaydataind_colHCD};
                for kk = 1:length(theofrag_HCD)
                    theofrag_HCD(kk).isomerind = jj;
                end
                theofrag_HCDsto{jj} = theofrag_HCD;                                     % stores HCD fragmentation spectra for given isomer
            end
            if any(displaydataind_colCID)
                theofrag_CID = temptheofrag{1,displaydataind_colCID};
                for kk = 1:length(theofrag_CID)
                    theofrag_CID(kk).isomerind = jj;
                end
                theofrag_CIDsto{jj} = theofrag_CID;
            end
            if any(displaydataind_colEThcD)
                theofrag_EThcD_cofrag = temptheofrag{1,displaydataind_colEThcD};
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
            if any(displaydataind_colETD)
                theofrag_ETD = temptheofrag{1,displaydataind_colETD};
                for kk = 1:length(theofrag_ETD)
                    theofrag_ETD(kk).isomerind = jj;
                end
                theofrag_ETDsto{jj} = theofrag_ETD;
            end
            if any(displaydataind_colETciD)
                theofrag_ETciD = temptheofrag{1,displaydataind_colETciD};
                for kk = 1:length(theofrag_ETciD)
                    theofrag_ETciD(kk).isomerind = jj;
                end
                theofrag_ETciDsto{jj} = theofrag_ETciD;
            end
        end


        %% CLEAN UP ISOMER LIST - DIAGNOSTIC ION BASED
        % Search for diagnostic ion, if found any, remove isomers that
        % don't contain this struct feature
        %         diagnosticmass_CID = [];
        %         diagnosticmasstyp_CID = [];
        %         diagnosticmass_HCD = [];
        %         diagnosticmasstyp_HCD = [];
        diagnosticmass_B = [];
        diagnosticmasstyp_B = [];
        diagnosticmasscount_B = zeros(numspecialfeatures,1);
        diagnosticmass_Y = [];
        diagnosticmasstyp_Y = [];
        diagnosticmasscount_Y = zeros(numspecialfeatures,1);
        %         diagnosticmass_Blacdinac = [];
        %         diagnosticmasstyp_Blacdinac = [];
        %         diagnosticmasscount_Blacdinac = zeros(numspecialfeatures,1);
        %         diagnosticmass_Ylacdinac = [];
        %         diagnosticmasstyp_Ylacdinac = [];
        %         diagnosticmasscount_Ylacdinac = zeros(numspecialfeatures,1);
        if any(structdetail_specialstruct(1,:))
            g_diag = g_ori;
            g_diag.struct = '{n{f}}';
            diag_1 = glypepMW(joinGlyPep(p,g_diag,m)) + 1.007825032;
            g_diag.struct = '{n{f}{n}}';
            diag_2 = glypepMW(joinGlyPep(p,g_diag,m)) + 1.007825032;
            diagnosticmass_Y = [diagnosticmass_Y;diag_1;diag_2];        % write glycan mass for core Fucose
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;1;1];
            diagnosticmasscount_Y(1) = 2;
        end
        if any(structdetail_specialstruct(2,:))
            g_diag = g_ori;
            g_diag.struct = '{n{n{h{n}}}}';
            diag_1 = glypepMW(joinGlyPep(p,g_diag,m)) + 1.007825032;
            g_diag.struct = '{n{f}{n{h{n}}}}';
            diag_2 = glypepMW(joinGlyPep(p,g_diag,m)) + 1.007825032;    % write glycan mass for bisecting
            diagnosticmass_Y = [diagnosticmass_Y;diag_1;diag_2];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;2;2];
            diagnosticmasscount_Y(2) = 1;
        end
        if any(structdetail_specialstruct(3,:))                         % LeX/a
            diag_1 = glyMW('{n{f}{h}}') + 1.007825032 - 18.0105647;
%            diag_2 = glyMW('{n{f}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;3];
            diagnosticmasscount_B(3) = 1;
            diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*1.007825032];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;3];
            diagnosticmasscount_Y(3) = 1;
            %             diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*(1.007825032 - 18.0105647);...
            %                 gpmw - diag_2 + 2*(1.007825032 - 18.0105647)];
            %             diagnosticmasstyp_Y = [diagnosticmasstyp_Y;3;3];
            %             diagnosticmasscount_Y(3) = 2;
        end
        if any(structdetail_specialstruct(4,:))                         % sLeX/a 
            diag_1 = glyMW('{n{f}{h{s}}}') + 1.007825032 - 18.0105647;
%            diag_2 = glyMW('{n{f}}') + 1.007825032 - 18.0105647;
%            diag_3 = glyMW('{n{h{s}}}') + 1.007825032 - 18.0105647;
%            diagnosticmass_B = [diagnosticmass_B;diag_1;diag_2;diag_3];
%            diagnosticmasstyp_B = [diagnosticmasstyp_B;4;4;4];
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;4];
            diagnosticmasscount_B(4) = 1;       % how many have to be found to say that the diagnostic ion exists
            diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*1.007825032];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;4];
            diagnosticmasscount_Y(4) = 1;
            %             diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*(1.007825032 - 18.0105647);...
            %                 gpmw - diag_2 + 2*(1.007825032 - 18.0105647);...
            %                 gpmw - diag_3 + 2*(1.007825032 - 18.0105647)];
            %             diagnosticmasstyp_Y = [diagnosticmasstyp_Y;4;4;4];
            %             diagnosticmasscount_Y(4) = 3;
        end
        if any(structdetail_specialstruct(5,:))                     % terminal sialic acid
            diag_1 = glyMW('{s}') + 1.007825032 - 18.0105647;
            diag_2 = glyMW('{n{h{s}}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1;diag_2];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;5;5];
            diagnosticmasscount_B(5) = 1;
            diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*1.007825032];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;5];
            diagnosticmasscount_Y(5) = 1;
            %             diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*(1.007825032 - 18.0105647);...
            %                 gpmw - diag_2 + 2*(1.007825032 - 18.0105647)];
            %             diagnosticmasstyp_Y = [diagnosticmasstyp_Y;5;5];
            %             diagnosticmasscount_Y(5) = 2;
        end
        if any(structdetail_specialstruct(6,:))                     % LacdiNAc
            diag_1 = glyMW('{n{n}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;6];
            diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*1.007825032];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;6];
            diagnosticmasscount_B(6) = 1;
            diagnosticmasscount_Y(6) = 1;
            %             tempdiagnosticmass_Blacdinac = [];
            %             tempdiagnosticmass_Ylacdinac = [];
            %             for jj = 1:length(isomers)
            %                 bondmap = getglycanbondmap(isomers{jj});
            %                 monosac = regexp(isomers{jj},'[hnsf]','match');
            %                 ndist = distance(strcmpi(monosac,'n'));
            %                 antennarootNind = find(ndist == 5);
            %                 for kk = 1:length(antennarootNind)
            %                     Nchildren = glytreetracker(bondmap,antennarootNind(kk),[],'down');
            %                     if ismember('n',monosac(Nchildren))
            %                         tempdiagnosticmass_Blacdinac = [tempdiagnosticmass_Blacdinac;...
            %                             glyMW(strjoin(monosac(Nchildren))) + 1.007825032 - 18.0105647];
            %                         tempdiagnosticmass_Ylacdinac = [tempdiagnosticmass_Ylacdinac;...
            %                             glyMW(isomers{jj}) - glyMW(strjoin(monosac(Nchildren))) + 1.007825032 + 18.0105647];
            %                     end
            %                 end
            %             end
            %             diagnosticmass_Blacdinac = [diagnosticmass_Blacdinac;tempdiagnosticmass_Blacdinac];
            %             diagnosticmasstyp_Blacdinac = [diagnosticmasstyp_Blacdinac;ones(size(tempdiagnosticmass_Blacdinac)) * 6];
            %             diagnosticmasscount_Blacdinac(6) = length(tempdiagnosticmass_Blacdinac);
            %             diagnosticmass_Ylacdinac = [diagnosticmass_Ylacdinac;tempdiagnosticmass_Ylacdinac];
            %             diagnosticmasstyp_Ylacdinac = [diagnosticmasstyp_Ylacdinac;ones(size(tempdiagnosticmass_Ylacdinac)) * 6];
            %             diagnosticmasscount_Ylacdinac(6) = length(tempdiagnosticmass_Ylacdinac);
        end
        if any(structdetail_specialstruct(7,:))                 % blood gp A
            diag_1 = glyMW('{n{h{f}{n}}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;7];
            diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*1.007825032];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;7];
            diagnosticmasscount_B(7) = 1;
            diagnosticmasscount_Y(7) = 1;
        end
        if any(structdetail_specialstruct(8,:))             % blood gp B
            diag_1 = glyMW('{h{f}{h}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;8];
            diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*1.007825032];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;8];
            diagnosticmasscount_B(8) = 1;
            diagnosticmasscount_Y(8) = 1;
        end
        if any(structdetail_specialstruct(9,:))             % blood gp O 
            diag_1 = glyMW('{h{f}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;9];
            diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*1.007825032];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;9];
            diagnosticmasscount_B(9) = 1;
            diagnosticmasscount_Y(9) = 1;
        end
        if any(structdetail_specialstruct(10,:))            % Lewis-Y/b
            diag_1 = glyMW('{n{f}{h{f}}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;10];
            diagnosticmasscount_B(10) = 1;
            diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*1.007825032];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;10];
            diagnosticmasscount_Y(10) = 1;
            %             diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*(1.007825032 - 18.0105647);...
            %                 gpmw - diag_2 + 2*(1.007825032 - 18.0105647);...
            %                 gpmw - diag_3 + 2*(1.007825032 - 18.0105647)];
            %             diagnosticmasstyp_Y = [diagnosticmasstyp_Y;10;10;10];
            %             diagnosticmasscount_Y(10) = 3;
        end
        if any(structdetail_specialstruct(11,:))        % Fucosylated LacdiNAc
            diag_1 = glyMW('{n{f}{n}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1;diag_2;diag_3];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;11];
            diagnosticmasscount_B(11) = 1;
            diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*1.007825032];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;11];
            diagnosticmasscount_Y(11) = 1;
            %             diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*(1.007825032 - 18.0105647);...
            %                 gpmw - diag_2 + 2*(1.007825032 - 18.0105647);...
            %                 gpmw - diag_3 + 2*(1.007825032 - 18.0105647)];
            %             diagnosticmasstyp_Y = [diagnosticmasstyp_Y;11;11;11];
            %             diagnosticmasscount_Y(11) = 3;
        end
        if any(structdetail_specialstruct(12,:))        % Sialylated LacdiNAc
            diag_1 = glyMW('{n{n{s}}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;12];
            diagnosticmasscount_B(12) = 1;
            diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*1.007825032];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;12];
            diagnosticmasscount_Y(12) = 1;
            %             diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*(1.007825032 - 18.0105647);...
            %                 gpmw - diag_2 + 2*(1.007825032 - 18.0105647);...
            %                 gpmw - diag_3 + 2*(1.007825032 - 18.0105647);...
            %                 gpmw - diag_4 + 2*(1.007825032 - 18.0105647)];
            %             diagnosticmasstyp_Y = [diagnosticmasstyp_Y;12;12;12;12];
            %             diagnosticmasscount_Y(12) = 4;
        end
        if any(structdetail_specialstruct(13,:))        % hybrid
            diag_1 = glyMW('{h{h{h}}}') + 1.007825032 - 18.0105647;
            diag_2 = glyMW('{h{h}}') + 1.007825032 - 18.0105647;
            diagnosticmass_B = [diagnosticmass_B;diag_1;diag_2];
            diagnosticmasstyp_B = [diagnosticmasstyp_B;13;13];
            diagnosticmasscount_B(13) = 1;
            diagnosticmass_Y = [diagnosticmass_Y;gpmw - diag_1 + 2*1.007825032];
            diagnosticmasstyp_Y = [diagnosticmasstyp_Y;13];
            diagnosticmasscount_Y(13) = 1;
        end
        if ~isempty(diagnosticmass_B)
            if any(methodind_CID)
                found_specialstruct_CID = false(numspecialfeatures,1);
                diagmatch_CID_B = calcithscore(spectrumCID,diagnosticmass_B,1,...
                    ms2tol(methodind_CID),ms2tolunit{methodind_CID},3,struct('maxlag',1,'selectpeak',[]));   % match theoretical fragments to experimental data.
                diagmatch_CID_Y = calcithscore(spectrumCID,diagnosticmass_Y,chgCID,...
                    ms2tol(methodind_CID),ms2tolunit{methodind_CID},3,struct('maxlag',1,'selectpeak',[]));
                %                 diagmatch_CID_Blacdinac = calcithscore(spectrumCID,diagnosticmass_Blacdinac,1,...
                %                     ms2tol(methodind_CID),ms2tolunit{methodind_CID},3,struct('maxlag',1,'selectpeak',[]));
                %                 diagmatch_CID_Ylacdinac = calcithscore(spectrumCID,diagnosticmass_Ylacdinac,1,...
                %                     ms2tol(methodind_CID),ms2tolunit{methodind_CID},3,struct('maxlag',1,'selectpeak',[]));
                matcheddiagiontyp_CID = [diagnosticmasstyp_B(logical(diagmatch_CID_B.ionmatchindex));...
                    diagnosticmasstyp_Y(logical(diagmatch_CID_Y.ionmatchindex))];
                %                 matcheddiagiontyp_CID_Blacdinac = logical(diagmatch_CID_Blacdinac.ionmatchindex);
                %                 matcheddiagiontyp_CID_Ylacdinac = logical(diagmatch_CID_Ylacdinac.ionmatchindex);
                for jj = 1:numspecialfeatures
                    if diagionlogicandor == 1   % both b and y ions
                        %                         if jj == 6
                        %                             featurefoundbydiagion = matcheddiagiontyp_CID_Blacdinac & matcheddiagiontyp_CID_Ylacdinac;
                        %                             if sum(featurefoundbydiagion) ==  diagnosticmasscount_Blacdinac(6) % found feature 1~9
                        %                                 found_specialstruct_CID(jj) = true;
                        %                             end
                        if any(structdetail_specialstruct(jj,:) == 1)
                            % mixed - isomer may or maynot has struct feature
                            featurefoundbydiagion_B = diagnosticmasstyp_B(logical(diagmatch_CID_B.ionmatchindex)) == jj;
                            featurefoundbydiagion_Y = diagnosticmasstyp_Y(logical(diagmatch_CID_Y.ionmatchindex)) == jj;
                            if (sum(featurefoundbydiagion_B) >= diagnosticmasscount_B(jj)) && ...
                                    (sum(featurefoundbydiagion_Y) >= diagnosticmasscount_Y(jj)) % found feature 1~9
                                found_specialstruct_CID(jj) = true;
                            end
                        end
                    elseif diagionlogicandor == 2   % either b or y ions
                        %                         if jj == 6
                        %                             featurefoundbydiagion = matcheddiagiontyp_CID_Blacdinac & matcheddiagiontyp_CID_Ylacdinac;
                        %                             if any(featurefoundbydiagion)  % found feature 1~9
                        %                                 found_specialstruct_CID(jj) = true;
                        %                             end
                        if any(structdetail_specialstruct(jj,:) == 1)
                            % mixed - isomer may or maynot has struct feature
                            featurefoundbydiagion = matcheddiagiontyp_CID == jj;
                            if any(featurefoundbydiagion)  % found feature 1~9
                                found_specialstruct_CID(jj) = true;
                            end
                        end
                    end
                end
            end
        else  % no feature to be found
            found_specialstruct_CID = true(numspecialfeatures,1);
        end
        if ~isempty(diagnosticmass_B)
            if any(methodind_HCD)
                found_specialstruct_HCD = false(numspecialfeatures,1);
                diagmatch_HCD_B = calcithscore(spectrumHCD,diagnosticmass_B,1,...
                    ms2tol(methodind_HCD),ms2tolunit{methodind_HCD},3,struct('maxlag',1,'selectpeak',[]));
                diagmatch_HCD_Y = calcithscore(spectrumHCD,diagnosticmass_Y,chgHCD,...
                    ms2tol(methodind_HCD),ms2tolunit{methodind_HCD},3,struct('maxlag',1,'selectpeak',[]));
                %                 diagmatch_HCD_Blacdinac = calcithscore(spectrumHCD,diagnosticmass_Blacdinac,1,...
                %                     ms2tol(methodind_HCD),ms2tolunit{methodind_HCD},3,struct('maxlag',1,'selectpeak',[]));
                %                 diagmatch_HCD_Ylacdinac = calcithscore(spectrumHCD,diagnosticmass_Ylacdinac,1,...
                %                     ms2tol(methodind_HCD),ms2tolunit{methodind_HCD},3,struct('maxlag',1,'selectpeak',[]));
                matcheddiagiontyp_HCD = [diagnosticmasstyp_B(logical(diagmatch_HCD_B.ionmatchindex));...
                    diagnosticmasstyp_Y(logical(diagmatch_HCD_Y.ionmatchindex))];
                %                 matcheddiagiontyp_HCD_Blacdinac = logical(diagmatch_HCD_Blacdinac.ionmatchindex);
                %                 matcheddiagiontyp_HCD_Ylacdinac = logical(diagmatch_HCD_Ylacdinac.ionmatchindex);
                for jj = 1:numspecialfeatures
                    if diagionlogicandor == 1   % AND
                        if any(structdetail_specialstruct(jj,:) == 1)
                            featurefoundbydiagion_B = diagnosticmasstyp_B(logical(diagmatch_HCD_B.ionmatchindex)) == jj;
                            featurefoundbydiagion_Y = diagnosticmasstyp_Y(logical(diagmatch_HCD_Y.ionmatchindex)) == jj;
                            if (sum(featurefoundbydiagion_B) >= diagnosticmasscount_B(jj)) && ...
                                    (sum(featurefoundbydiagion_Y) >= diagnosticmasscount_Y(jj)) % found feature 1~9
                                found_specialstruct_HCD(jj) = true;
                            end
                        end
                    elseif diagionlogicandor == 2   % OR
                        if any(structdetail_specialstruct(jj,:) == 1)
                            featurefoundbydiagion = matcheddiagiontyp_HCD == jj;
                            if any(featurefoundbydiagion)  % found feature 1~9
                                found_specialstruct_HCD(jj) = true;
                            end
                        end
                    end
                end
            end
        else
            found_specialstruct_HCD = true(numspecialfeatures,1);
        end
        if ~isempty(diagnosticmass_B)
            if any(methodind_EThcD)
                found_specialstruct_EThcD = false(numspecialfeatures,1);
                diagmatch_EThcD_B = calcithscore(spectrumEThcD,diagnosticmass_B,1,...
                    ms2tol(methodind_EThcD),ms2tolunit{methodind_EThcD},3,struct('maxlag',1,'selectpeak',[]));
                diagmatch_EThcD_Y = calcithscore(spectrumEThcD,diagnosticmass_Y,chgEThcD,...
                    ms2tol(methodind_EThcD),ms2tolunit{methodind_EThcD},3,struct('maxlag',1,'selectpeak',[]));
                %                 diagmatch_EThcD_Blacdinac = calcithscore(spectrumEThcD,diagnosticmass_Blacdinac,1,...
                %                     ms2tol(methodind_EThcD),ms2tolunit{methodind_EThcD},3,struct('maxlag',1,'selectpeak',[]));
                %                 diagmatch_EThcD_Ylacdinac = calcithscore(spectrumEThcD,diagnosticmass_Ylacdinac,1,...
                %                     ms2tol(methodind_EThcD),ms2tolunit{methodind_EThcD},3,struct('maxlag',1,'selectpeak',[]));
                matcheddiagiontyp_EThcD = [diagnosticmasstyp_B(logical(diagmatch_EThcD_B.ionmatchindex));...
                    diagnosticmasstyp_Y(logical(diagmatch_EThcD_Y.ionmatchindex))];
                %                 matcheddiagiontyp_EThcD_Blacdinac = logical(diagmatch_EThcD_Blacdinac.ionmatchindex);
                %                 matcheddiagiontyp_EThcD_Ylacdinac = logical(diagmatch_EThcD_Ylacdinac.ionmatchindex);
                for jj = 1:numspecialfeatures
                    if diagionlogicandor == 1
                        %                         if jj == 6
                        %                             featurefoundbydiagion = matcheddiagiontyp_EThcD_Blacdinac & matcheddiagiontyp_EThcD_Ylacdinac;
                        %                             if sum(featurefoundbydiagion) ==  diagnosticmasscount_Blacdinac(6) % found feature 1~9
                        %                                 found_specialstruct_EThcD(jj) = true;
                        %                             end
                        if any(structdetail_specialstruct(jj,:) == 1)
                            featurefoundbydiagion_B = diagnosticmasstyp_B(logical(diagmatch_EThcD_B.ionmatchindex)) == jj;
                            featurefoundbydiagion_Y = diagnosticmasstyp_Y(logical(diagmatch_EThcD_Y.ionmatchindex)) == jj;
                            if (sum(featurefoundbydiagion_B) >= diagnosticmasscount_B(jj)) && ...
                                    (sum(featurefoundbydiagion_Y) >= diagnosticmasscount_Y(jj)) % found feature 1~9
                                found_specialstruct_EThcD(jj) = true;
                            end
                        end
                    elseif diagionlogicandor == 2
                        if any(structdetail_specialstruct(jj,:) == 1)
                            % mixed - isomer may or maynot has struct feature
                            featurefoundbydiagion = matcheddiagiontyp_EThcD == jj;
                            if any(featurefoundbydiagion)  % found feature 1~9
                                found_specialstruct_EThcD(jj) = true;
                            end
                        end
                    end
                end
            end
        else
            found_specialstruct_EThcD = true(numspecialfeatures,1);
        end
        if ~isempty(diagnosticmass_B)
            if any(methodind_ETciD)
                found_specialstruct_ETciD = false(numspecialfeatures,1);
                diagmatch_ETciD_B = calcithscore(spectrumETciD,diagnosticmass_B,1,...
                    ms2tol(methodind_ETciD),ms2tolunit{methodind_ETciD},3,struct('maxlag',1,'selectpeak',[]));
                diagmatch_ETciD_Y = calcithscore(spectrumETciD,diagnosticmass_Y,chgETciD,...
                    ms2tol(methodind_ETciD),ms2tolunit{methodind_ETciD},3,struct('maxlag',1,'selectpeak',[]));
                matcheddiagiontyp_ETciD = [diagnosticmasstyp_B(logical(diagmatch_ETciD_B.ionmatchindex));...
                    diagnosticmasstyp_Y(logical(diagmatch_ETciD_Y.ionmatchindex))];
                for jj = 1:numspecialfeatures
                    if diagionlogicandor == 1
                        if any(structdetail_specialstruct(jj,:) == 1)
                            % mixed - isomer may or maynot has struct feature
                            featurefoundbydiagion_B = diagnosticmasstyp_B(logical(diagmatch_ETciD_B.ionmatchindex)) == jj;
                            featurefoundbydiagion_Y = diagnosticmasstyp_Y(logical(diagmatch_ETciD_Y.ionmatchindex)) == jj;
                            if (sum(featurefoundbydiagion_B) >= diagnosticmasscount_B(jj)) && ...
                                    (sum(featurefoundbydiagion_Y) >= diagnosticmasscount_Y(jj)) % found feature 1~9
                                found_specialstruct_ETciD(jj) = true;
                            end
                        end
                    elseif diagionlogicandor == 2
                        if any(structdetail_specialstruct(jj,:) == 1)
                            % mixed - isomer may or maynot has struct feature
                            featurefoundbydiagion = matcheddiagiontyp_ETciD == jj;
                            if any(featurefoundbydiagion)  % found feature 1~9
                                found_specialstruct_ETciD(jj) = true;
                            end
                        end
                    end
                end
            end
        else  % no feature to be found
            found_specialstruct_ETciD = true(numspecialfeatures,1);
        end
        found_specialstruct_ETD = true(numspecialfeatures,1);   % not used since we are not setup for ETD
        isomerkeepind = false(size(isomers));
        isomerhasstructfeatures = any(any(structdetail_specialstruct));
        if isomerhasstructfeatures  % any isomer has structfeatures
            for jj = 1:length(isomers)
                if any(structdetail_specialstruct(:,jj))  % this isomer has some structfeatures
                    found_specialstruct = zeros(numspecialfeatures,1);
                    if any(displaydataind_colCID)
                        found_specialstruct = found_specialstruct + double(found_specialstruct_CID);
                    end
                    if any(displaydataind_colHCD)
                        found_specialstruct = found_specialstruct + double(found_specialstruct_HCD);
                    end
                    if any(displaydataind_colEThcD)
                        found_specialstruct = found_specialstruct + double(found_specialstruct_EThcD);
                    end
                    if any(displaydataind_colETciD)
                        found_specialstruct = found_specialstruct + double(found_specialstruct_ETciD);
                    end
                    keepthisisomer = ~xor(structdetail_specialstruct(:,jj),...
                        structdetail_specialstruct(:,jj) & (found_specialstruct >= nummindiagionfragmode));  % feature found
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
            if ~any(isomerkeepind)      % if no isomers are found then keep all
                isomerkeepind = true(size(isomers));
            end
        else  % no isomer has structfeatures, then keep all isomers
            isomerkeepind = true(size(isomers));
        end
        isomers = isomers(isomerkeepind);
        if isempty(isomers)
            outputtoresult = false;
            ct=ct+1;
            ii, ct
        else
            theofrag_HCDsto = theofrag_HCDsto(isomerkeepind);
            theofrag_EThcD_cofragsto = theofrag_EThcD_cofragsto(isomerkeepind);
            theofrag_EThcD_glyonlysto = theofrag_EThcD_glyonlysto(isomerkeepind);
            theofrag_CIDsto = theofrag_CIDsto(isomerkeepind);
            theofrag_ETciDsto = theofrag_ETciDsto(isomerkeepind);
            isomerProteinIDsto = isomerProteinIDsto(isomerkeepind);                 % This just stores the IDs for protein and glycans for later use
            for jj = 1:length(isomers)
                [theofrag_HCDsto{jj}.isomerind] = deal(jj);
                [theofrag_EThcD_cofragsto{jj}.isomerind] = deal(jj);
                [theofrag_EThcD_glyonlysto{jj}.isomerind] = deal(jj);
                [theofrag_CIDsto{jj}.isomerind] = deal(jj);
                [theofrag_ETciDsto{jj}.isomerind] = deal(jj);
            end

            for jj = 1:length(isomers)
                isomerfrags_HCD = [isomerfrags_HCD,theofrag_HCDsto{jj}];
                isomerfrags_EThcD_glyonly = [isomerfrags_EThcD_glyonly,theofrag_EThcD_glyonlysto{jj}];
                isomerfrags_CID = [isomerfrags_CID,theofrag_CIDsto{jj}];
                isomerfrags_ETciD = [isomerfrags_ETciD,theofrag_ETciDsto{jj}];
            end

            %% MATCH EXPT SPECTRUM
            if any(displaydataind_colHCD)
                isomerfrags_HCDmz = [];
                isomerfrags_HCDmzisomerind = [];
                isomerfrags_HCDchg = [];
                isomerfrags_HCDfragind = [];
                isomerfrags_HCD_originalisomerind = [isomerfrags_HCD.isomerind];
                for jj = 1:chgHCD  % get the biglist
                    tempisomerfrags_HCDmz = ([isomerfrags_HCD.mz] - 1.007825032)/jj + 1.007825032;
                    isomerfrags_HCDmz = [isomerfrags_HCDmz;tempisomerfrags_HCDmz(:)];
                    isomerfrags_HCDmzisomerind = [isomerfrags_HCDmzisomerind;[isomerfrags_HCD.isomerind]];
                    isomerfrags_HCDchg = [isomerfrags_HCDchg;ones(length(isomerfrags_HCD),1) * jj];
                    isomerfrags_HCDfragind = [isomerfrags_HCDfragind;(1:length(isomerfrags_HCD))'];
                end
                [~,~,HCDisomerfrags_uniqueness] = unique(isomerfrags_HCDmz,'stable');
                %% SPECTRUM MATCHING - HCD - GLYONLY
                match_HCD = calcithscore(spectrumHCD,[isomerfrags_HCD.mz],chgHCD,...
                    ms2tol(methodind_HCD),ms2tolunit{methodind_HCD},2,...
                    struct('maxlag',scoreoptions.maxlag,'selectpeak',[]));
                medianpeakHCD = median(spectrumHCD(:,2));
                matchedpeakind_HCD = match_HCD.peakmatchindex;
                matchedpeakind_HCD_type = zeros(size(matchedpeakind_HCD,1),1);  % temp variable, store info for potential use
                HCD_uniisomers_numerical = zeros(1,length(isomers));
                HCD_sharedisomers_numerical = zeros(1,length(isomers));
                HCD_uniisomers_quant = zeros(1,length(isomers));
                HCD_sharedisomers_quant = zeros(1,length(isomers));
                allmatchedpeakind_HCD = [];
                allmatchedpeakind_HCD_whichpeak = [];
                for jj = 1:size(matchedpeakind_HCD,1)
                    if ~isempty(matchedpeakind_HCD{jj})
                        allmatchedpeakind_HCD = [allmatchedpeakind_HCD;matchedpeakind_HCD{jj}];
                        allmatchedpeakind_HCD_whichpeak = [allmatchedpeakind_HCD_whichpeak;repmat(jj,size(matchedpeakind_HCD{jj},1),1)];
                    end
                end
                [allmatchedpeakind_HCD,~,allmatchedpeakind_HCD_uniind] = unique(allmatchedpeakind_HCD,'rows');
                tempHCDimportantstructfeaturesfound_isomerid = cell(size(allmatchedpeakind_HCD,1),1);
                for jj = 1:size(allmatchedpeakind_HCD,1)  % fragind, chg
                    tempHCDmatchedunifragind = allmatchedpeakind_HCD(jj,1);  % matched frag ind
                    tempHCDmatchedunifragchg = allmatchedpeakind_HCD(jj,2);  % maeched frag chg
                    tempHCDmatchedpeakind = unique(allmatchedpeakind_HCD_whichpeak(allmatchedpeakind_HCD_uniind == jj));
                    tempHCDmatchedunifragind_inbig = isomerfrags_HCDfragind == tempHCDmatchedunifragind & ...
                        isomerfrags_HCDchg == tempHCDmatchedunifragchg;  % find its absolute position in biglist
                    tempHCDmatchedfraguniqueind = HCDisomerfrags_uniqueness(tempHCDmatchedunifragind_inbig);
                    matchedfragind_HCD_fororiginal = isomerfrags_HCDfragind(ismember(HCDisomerfrags_uniqueness,tempHCDmatchedfraguniqueind));
                    HCDmatchedorifrag  = isomerfrags_HCD(matchedfragind_HCD_fororiginal);
                    HCDmatchedorifragchg = isomerfrags_HCDchg(ismember(HCDisomerfrags_uniqueness,tempHCDmatchedfraguniqueind));
                    keepthisHCDfrag = true(size(HCDmatchedorifrag));
                    sigionquant_HCD = sum(log2(spectrumHCD(tempHCDmatchedpeakind,2)/medianpeakHCD));
                    for kk = 1:length(HCDmatchedorifrag)
                        temptype = HCDmatchedorifrag(kk).type;
                        if strcmpi(temptype(1:2),'-b')  % Bion
                            if HCDmatchedorifragchg(kk) > 1
                                keepthisHCDfrag(kk) = false;
                            end
                        elseif strcmpi(temptype(1:2),'-y')  % Bion
                            if HCDmatchedorifragchg(kk) > 3
                                keepthisHCDfrag(kk) = false;
                            end
                        end
                    end
                    matchedfragind_HCD_fororiginal = matchedfragind_HCD_fororiginal(keepthisHCDfrag);
                    [tempmatchedisomerind,~,tempmatchedisomerindind2] = ...
                        unique(isomerfrags_HCD_originalisomerind(matchedfragind_HCD_fororiginal));
                    if length(tempmatchedisomerind) == 1  % unique
                        matchedpeakind_HCD_type(jj) = 1;
                        HCD_uniisomers_numerical(tempmatchedisomerind) = HCD_uniisomers_numerical(tempmatchedisomerind) + ...
                            length(tempmatchedisomerindind2);
                        HCD_uniisomers_quant(tempmatchedisomerind) = HCD_uniisomers_quant(tempmatchedisomerind) + ...
                            sigionquant_HCD;
                        tempHCDimportantstructfeaturesfound_isomerid{jj} = tempmatchedisomerind;
                    elseif length(tempmatchedisomerind) == length(isomers)  % common
                        matchedpeakind_HCD_type(jj) = -1;
                    elseif ~isempty(tempmatchedisomerind)  % shared
                        matchedpeakind_HCD_type(jj) = 2;
                        tempHCDisomerind = [];
                        for kk = 1:length(tempmatchedisomerind)
                            numisomerfound = sum(tempmatchedisomerindind2 == tempmatchedisomerindind2(kk));
                            HCD_sharedisomers_numerical(tempmatchedisomerind(kk)) = ...
                                HCD_sharedisomers_numerical(tempmatchedisomerind(kk)) + numisomerfound;
                            HCD_sharedisomers_quant(tempmatchedisomerind(kk)) = HCD_sharedisomers_quant(tempmatchedisomerind(kk)) + ...
                                sigionquant_HCD;
                            tempHCDisomerind = [tempHCDisomerind,tempmatchedisomerind(kk)];
                        end
                    end
                end
            end
            if any(displaydataind_colEThcD)
                isomerfrags_EThcD_glyonlymz = [];
                isomerfrags_EThcD_glyonlymzisomerind = [];
                isomerfrags_EThcD_glyonlychg = [];
                isomerfrags_EThcD_glyonlyfragind = [];
                isomerfrags_EThcD_originalisomerind = [isomerfrags_EThcD_glyonly.isomerind];
                for jj = 1:chgEThcD  % get the biglist
                    tempisomerfrags_EThcD_glyonlymz = ([isomerfrags_EThcD_glyonly.mz] - 1.007825032)/jj + 1.007825032;
                    isomerfrags_EThcD_glyonlymz = [isomerfrags_EThcD_glyonlymz;tempisomerfrags_EThcD_glyonlymz(:)];
                    isomerfrags_EThcD_glyonlymzisomerind = [isomerfrags_EThcD_glyonlymzisomerind;[isomerfrags_EThcD_glyonly.isomerind]];
                    isomerfrags_EThcD_glyonlychg = [isomerfrags_EThcD_glyonlychg;ones(length(isomerfrags_EThcD_glyonly),1) * jj];
                    isomerfrags_EThcD_glyonlyfragind = [isomerfrags_EThcD_glyonlyfragind;(1:length(isomerfrags_EThcD_glyonly))'];
                end
                [~,~,EThcDisomerfrags_uniqueness] = unique(isomerfrags_EThcD_glyonlymz,'stable');
                match_EThcD = calcithscore(spectrumEThcD,[isomerfrags_EThcD_glyonly.mz],chgEThcD,...
                    ms2tol(methodind_EThcD),ms2tolunit{methodind_EThcD},2,...
                    struct('maxlag',scoreoptions.maxlag,'selectpeak',[]));
                medianpeakEThcD = median(spectrumEThcD(:,2));
                matchedpeakind_EThcD = match_EThcD.peakmatchindex;
                matchedpeakind_EThcD_type = zeros(size(matchedpeakind_EThcD,1),1);  % temp variable, store info for potential use
                EThcD_uniisomers_numerical = zeros(1,length(isomers));
                EThcD_sharedisomers_numerical = zeros(1,length(isomers));
                EThcD_uniisomers_quant = zeros(1,length(isomers));
                EThcD_sharedisomers_quant = zeros(1,length(isomers));
                allmatchedpeakind_EThcD = [];
                allmatchedpeakind_EThcD_whichpeak = [];
                for jj = 1:size(matchedpeakind_EThcD,1)
                    if ~isempty(matchedpeakind_EThcD{jj})
                        allmatchedpeakind_EThcD = [allmatchedpeakind_EThcD;matchedpeakind_EThcD{jj}];
                        allmatchedpeakind_EThcD_whichpeak = [allmatchedpeakind_EThcD_whichpeak;repmat(jj,size(matchedpeakind_EThcD{jj},1),1)];
                    end
                end
                [allmatchedpeakind_EThcD,~,allmatchedpeakind_EThcD_uniind] = unique(allmatchedpeakind_EThcD,'rows');
                for jj = 1:size(allmatchedpeakind_EThcD,1)  % fragind, chg
                    tempEThcDmatchedunifragind = allmatchedpeakind_EThcD(jj,1);  % matched frag ind
                    tempEThcDmatchedunifragchg = allmatchedpeakind_EThcD(jj,2);  % matched frag chg
                    tempEThcDmatchedpeakind = unique(allmatchedpeakind_EThcD_whichpeak(allmatchedpeakind_EThcD_uniind == jj));  % which peak got matched
                    tempEThcDmatchedunifragind_inbig = isomerfrags_EThcD_glyonlyfragind == tempEThcDmatchedunifragind & ...
                        isomerfrags_EThcD_glyonlychg == tempEThcDmatchedunifragchg;  % find its absolute position in biglist
                    tempEThcDmatchedfraguniqueind = EThcDisomerfrags_uniqueness(tempEThcDmatchedunifragind_inbig);
                    matchedfragind_EThcD_fororiginal = isomerfrags_EThcD_glyonlyfragind(ismember(EThcDisomerfrags_uniqueness,tempEThcDmatchedfraguniqueind));
                    EThcDmatchedorifrag  = isomerfrags_EThcD_glyonly(matchedfragind_EThcD_fororiginal);
                    EThcDmatchedorifragchg = isomerfrags_EThcD_glyonlychg(ismember(EThcDisomerfrags_uniqueness,tempEThcDmatchedfraguniqueind));
                    keepthisEThcDfrag = true(size(EThcDmatchedorifrag));
                    sigionquant_EThcD = sum(log2(spectrumEThcD(tempEThcDmatchedpeakind,2)/medianpeakEThcD));
                    for kk = 1:length(EThcDmatchedorifrag)
                        temptype = EThcDmatchedorifrag(kk).type;
                        if strcmpi(temptype(1:2),'-b')  % Bion
                            if EThcDmatchedorifragchg(kk) > 1
                                keepthisEThcDfrag(kk) = false;
                            end
                        end
                    end
                    matchedfragind_EThcD_fororiginal = matchedfragind_EThcD_fororiginal(keepthisEThcDfrag);
                    [tempmatchedisomerind,~,tempmatchedisomerindind2] = ...
                        unique(isomerfrags_EThcD_originalisomerind(matchedfragind_EThcD_fororiginal));
                    if length(tempmatchedisomerind) == 1  % unique
                        matchedpeakind_EThcD_type(jj) = 1;
                        EThcD_uniisomers_numerical(tempmatchedisomerind) = EThcD_uniisomers_numerical(tempmatchedisomerind) + ...
                            length(tempmatchedisomerindind2);
                        EThcD_uniisomers_quant(tempmatchedisomerind) = EThcD_uniisomers_quant(tempmatchedisomerind) + ...
                            sigionquant_EThcD;
                    elseif length(tempmatchedisomerind) == length(isomers)  % common
                        matchedpeakind_EThcD_type(jj) = -1;
                    elseif ~isempty(tempmatchedisomerind)  % shared
                        matchedpeakind_EThcD_type(jj) = 2;
                        tempEThcDisomerind = [];
                        for kk = 1:length(tempmatchedisomerind)
                            numisomerfound = sum(tempmatchedisomerindind2 == tempmatchedisomerindind2(kk));
                            EThcD_sharedisomers_numerical(tempmatchedisomerind(kk)) = ...
                                EThcD_sharedisomers_numerical(tempmatchedisomerind(kk)) + numisomerfound;
                            EThcD_sharedisomers_quant(tempmatchedisomerind(kk)) = EThcD_sharedisomers_quant(tempmatchedisomerind(kk)) + ...
                                sigionquant_EThcD;
                            tempEThcDisomerind = [tempEThcDisomerind,tempmatchedisomerind(kk)];
                        end
                    end
                end
            end
            if any(displaydataind_colCID)
                isomerfrags_CIDmz = [];
                isomerfrags_CIDmzisomerind = [];
                isomerfrags_CIDchg = [];
                isomerfrags_CIDfragind = [];
                isomerfrags_CID_originalisomerind = [isomerfrags_CID.isomerind];
                for jj = 1:chgCID  % get the biglist
                    tempisomerfrags_CIDmz = ([isomerfrags_CID.mz] - 1.007825032)/jj + 1.007825032;
                    isomerfrags_CIDmz = [isomerfrags_CIDmz;tempisomerfrags_CIDmz(:)];
                    isomerfrags_CIDmzisomerind = [isomerfrags_CIDmzisomerind;[isomerfrags_CID.isomerind]];
                    isomerfrags_CIDchg = [isomerfrags_CIDchg;ones(length(isomerfrags_CID),1) * jj];
                    isomerfrags_CIDfragind = [isomerfrags_CIDfragind;(1:length(isomerfrags_CID))'];
                end
                [~,~,CIDisomerfrags_uniqueness] = unique(isomerfrags_CIDmz,'stable');
                %% SPECTRUM MATCHING - CID - GLYONLY
                match_CID = calcithscore(spectrumCID,[isomerfrags_CID.mz],chgCID,...
                    ms2tol(methodind_CID),ms2tolunit{methodind_CID},2,...
                    struct('maxlag',scoreoptions.maxlag,'selectpeak',[]));
                medianpeakCID = median(spectrumCID(:,2));
                matchedpeakind_CID = match_CID.peakmatchindex;
                matchedpeakind_CID_type = zeros(size(matchedpeakind_CID,1),1);  % temp variable, store info for potential use
                CID_uniisomers_numerical = zeros(1,length(isomers));
                CID_sharedisomers_numerical = zeros(1,length(isomers));
                CID_uniisomers_quant = zeros(1,length(isomers));
                CID_sharedisomers_quant = zeros(1,length(isomers));
                allmatchedpeakind_CID = [];
                allmatchedpeakind_CID_whichpeak = [];
                for jj = 1:size(matchedpeakind_CID,1)
                    if ~isempty(matchedpeakind_CID{jj})
                        allmatchedpeakind_CID = [allmatchedpeakind_CID;matchedpeakind_CID{jj}];
                        allmatchedpeakind_CID_whichpeak = [allmatchedpeakind_CID_whichpeak;repmat(jj,size(matchedpeakind_CID{jj},1),1)];
                    end
                end
                [allmatchedpeakind_CID,~,allmatchedpeakind_CID_uniind] = unique(allmatchedpeakind_CID,'rows');
                tempCIDimportantstructfeaturesfound_isomerid = cell(size(allmatchedpeakind_CID,1),1);
                for jj = 1:size(allmatchedpeakind_CID,1)  % fragind, chg
                    tempCIDmatchedunifragind = allmatchedpeakind_CID(jj,1);  % matched frag ind
                    tempCIDmatchedunifragchg = allmatchedpeakind_CID(jj,2);  % maeched frag chg
                    tempCIDmatchedpeakind = unique(allmatchedpeakind_CID_whichpeak(allmatchedpeakind_CID_uniind == jj));
                    tempCIDmatchedunifragind_inbig = isomerfrags_CIDfragind == tempCIDmatchedunifragind & ...
                        isomerfrags_CIDchg == tempCIDmatchedunifragchg;  % find its absolute position in biglist
                    tempCIDmatchedfraguniqueind = CIDisomerfrags_uniqueness(tempCIDmatchedunifragind_inbig);
                    matchedfragind_CID_fororiginal = isomerfrags_CIDfragind(ismember(CIDisomerfrags_uniqueness,tempCIDmatchedfraguniqueind));
                    CIDmatchedorifrag  = isomerfrags_CID(matchedfragind_CID_fororiginal);
                    CIDmatchedorifragchg = isomerfrags_CIDchg(ismember(CIDisomerfrags_uniqueness,tempCIDmatchedfraguniqueind));
                    keepthisCIDfrag = true(size(CIDmatchedorifrag));
                    sigionquant_CID = sum(log2(spectrumCID(tempCIDmatchedpeakind,2)/medianpeakCID));
                    for kk = 1:length(CIDmatchedorifrag)
                        temptype = CIDmatchedorifrag(kk).type;
                        if strcmpi(temptype(1:2),'-b')  % Bion
                            if CIDmatchedorifragchg(kk) > 1
                                keepthisCIDfrag(kk) = false;
                            end
                        elseif strcmpi(temptype(1:2),'-y')  % Bion
                            if CIDmatchedorifragchg(kk) > 3
                                keepthisCIDfrag(kk) = false;
                            end
                        end
                    end
                    matchedfragind_CID_fororiginal = matchedfragind_CID_fororiginal(keepthisCIDfrag);
                    [tempmatchedisomerind,~,tempmatchedisomerindind2] = ...
                        unique(isomerfrags_CID_originalisomerind(matchedfragind_CID_fororiginal));
                    if length(tempmatchedisomerind) == 1  % unique
                        matchedpeakind_CID_type(jj) = 1;
                        CID_uniisomers_numerical(tempmatchedisomerind) = CID_uniisomers_numerical(tempmatchedisomerind) + ...
                            length(tempmatchedisomerindind2);
                        CID_uniisomers_quant(tempmatchedisomerind) = CID_uniisomers_quant(tempmatchedisomerind) + ...
                            sigionquant_CID;
                        tempCIDimportantstructfeaturesfound_isomerid{jj} = tempmatchedisomerind;
                    elseif length(tempmatchedisomerind) == length(isomers)  % common
                        matchedpeakind_CID_type(jj) = -1;
                    elseif ~isempty(tempmatchedisomerind)  % shared
                        matchedpeakind_CID_type(jj) = 2;
                        tempCIDisomerind = [];
                        for kk = 1:length(tempmatchedisomerind)
                            numisomerfound = sum(tempmatchedisomerindind2 == tempmatchedisomerindind2(kk));
                            CID_sharedisomers_numerical(tempmatchedisomerind(kk)) = ...
                                CID_sharedisomers_numerical(tempmatchedisomerind(kk)) + numisomerfound;
                            CID_sharedisomers_quant(tempmatchedisomerind(kk)) = CID_sharedisomers_quant(tempmatchedisomerind(kk)) + ...
                                sigionquant_CID;
                            tempCIDisomerind = [tempCIDisomerind,tempmatchedisomerind(kk)];
                        end
                    end
                end
            end
            if any(displaydataind_colETciD)
                isomerfrags_ETciDmz = [];
                isomerfrags_ETciDmzisomerind = [];
                isomerfrags_ETciDchg = [];
                isomerfrags_ETciDfragind = [];
                isomerfrags_ETciD_originalisomerind = [isomerfrags_ETciD.isomerind];
                for jj = 1:chgETciD  % get the biglist
                    tempisomerfrags_ETciDmz = ([isomerfrags_ETciD.mz] - 1.007825032)/jj + 1.007825032;
                    isomerfrags_ETciDmz = [isomerfrags_ETciDmz;tempisomerfrags_ETciDmz(:)];
                    isomerfrags_ETciDmzisomerind = [isomerfrags_ETciDmzisomerind;[isomerfrags_ETciD.isomerind]];
                    isomerfrags_ETciDchg = [isomerfrags_ETciDchg;ones(length(isomerfrags_ETciD),1) * jj];
                    isomerfrags_ETciDfragind = [isomerfrags_ETciDfragind;(1:length(isomerfrags_ETciD))'];
                end
                [~,~,ETciDisomerfrags_uniqueness] = unique(isomerfrags_ETciDmz,'stable');
                %% SPECTRUM MATCHING - ETciD
                match_ETciD = calcithscore(spectrumETciD,[isomerfrags_ETciD.mz],chgETciD,...
                    ms2tol(methodind_ETciD),ms2tolunit{methodind_ETciD},2,...
                    struct('maxlag',scoreoptions.maxlag,'selectpeak',[]));
                medianpeakETciD = median(spectrumETciD(:,2));
                matchedpeakind_ETciD = match_ETciD.peakmatchindex;
                matchedpeakind_ETciD_type = zeros(size(matchedpeakind_ETciD,1),1);  % temp variable, store info for potential use
                ETciD_uniisomers_numerical = zeros(1,length(isomers));
                ETciD_sharedisomers_numerical = zeros(1,length(isomers));
                ETciD_uniisomers_quant = zeros(1,length(isomers));
                ETciD_sharedisomers_quant = zeros(1,length(isomers));
                allmatchedpeakind_ETciD = [];
                allmatchedpeakind_ETciD_whichpeak = [];
                for jj = 1:size(matchedpeakind_ETciD,1)
                    if ~isempty(matchedpeakind_ETciD{jj})
                        allmatchedpeakind_ETciD = [allmatchedpeakind_ETciD;matchedpeakind_ETciD{jj}];
                        allmatchedpeakind_ETciD_whichpeak = [allmatchedpeakind_ETciD_whichpeak;repmat(jj,size(matchedpeakind_ETciD{jj},1),1)];
                    end
                end
                [allmatchedpeakind_ETciD,~,allmatchedpeakind_ETciD_uniind] = unique(allmatchedpeakind_ETciD,'rows');
                tempETciDimportantstructfeaturesfound_isomerid = cell(size(allmatchedpeakind_ETciD,1),1);
                for jj = 1:size(allmatchedpeakind_ETciD,1)  % fragind, chg
                    tempETciDmatchedunifragind = allmatchedpeakind_ETciD(jj,1);  % matched frag ind
                    tempETciDmatchedunifragchg = allmatchedpeakind_ETciD(jj,2);  % maeched frag chg
                    tempETciDmatchedpeakind = unique(allmatchedpeakind_ETciD_whichpeak(allmatchedpeakind_ETciD_uniind == jj));
                    tempETciDmatchedunifragind_inbig = isomerfrags_ETciDfragind == tempETciDmatchedunifragind & ...
                        isomerfrags_ETciDchg == tempETciDmatchedunifragchg;  % find its absolute position in biglist
                    tempETciDmatchedfraguniqueind = ETciDisomerfrags_uniqueness(tempETciDmatchedunifragind_inbig);
                    matchedfragind_ETciD_fororiginal = isomerfrags_ETciDfragind(ismember(ETciDisomerfrags_uniqueness,tempETciDmatchedfraguniqueind));
                    ETciDmatchedorifrag  = isomerfrags_ETciD(matchedfragind_ETciD_fororiginal);
                    ETciDmatchedorifragchg = isomerfrags_ETciDchg(ismember(ETciDisomerfrags_uniqueness,tempETciDmatchedfraguniqueind));
                    keepthisETciDfrag = true(size(ETciDmatchedorifrag));
                    sigionquant_ETciD = sum(log2(spectrumETciD(tempETciDmatchedpeakind,2)/medianpeakETciD));
                    for kk = 1:length(ETciDmatchedorifrag)
                        temptype = ETciDmatchedorifrag(kk).type;
                        if strcmpi(temptype(1:2),'-b')  % Bion
                            if ETciDmatchedorifragchg(kk) > 1
                                keepthisETciDfrag(kk) = false;
                            end
                        elseif strcmpi(temptype(1:2),'-y')  % Bion
                            if ETciDmatchedorifragchg(kk) > 3
                                keepthisETciDfrag(kk) = false;
                            end
                        end
                    end
                    matchedfragind_ETciD_fororiginal = matchedfragind_ETciD_fororiginal(keepthisETciDfrag);
                    [tempmatchedisomerind,~,tempmatchedisomerindind2] = ...
                        unique(isomerfrags_ETciD_originalisomerind(matchedfragind_ETciD_fororiginal));
                    if length(tempmatchedisomerind) == 1  % unique
                        matchedpeakind_ETciD_type(jj) = 1;
                        ETciD_uniisomers_numerical(tempmatchedisomerind) = ETciD_uniisomers_numerical(tempmatchedisomerind) + ...
                            length(tempmatchedisomerindind2);
                        ETciD_uniisomers_quant(tempmatchedisomerind) = ETciD_uniisomers_quant(tempmatchedisomerind) + ...
                            sigionquant_ETciD;
                        tempETciDimportantstructfeaturesfound_isomerid{jj} = tempmatchedisomerind;
                    elseif length(tempmatchedisomerind) == length(isomers)  % common
                        matchedpeakind_ETciD_type(jj) = -1;
                    elseif ~isempty(tempmatchedisomerind)  % shared
                        matchedpeakind_ETciD_type(jj) = 2;
                        tempETciDisomerind = [];
                        for kk = 1:length(tempmatchedisomerind)
                            numisomerfound = sum(tempmatchedisomerindind2 == tempmatchedisomerindind2(kk));
                            ETciD_sharedisomers_numerical(tempmatchedisomerind(kk)) = ...
                                ETciD_sharedisomers_numerical(tempmatchedisomerind(kk)) + numisomerfound;
                            ETciD_sharedisomers_quant(tempmatchedisomerind(kk)) = ETciD_sharedisomers_quant(tempmatchedisomerind(kk)) + ...
                                sigionquant_ETciD;
                            tempETciDisomerind = [tempETciDisomerind,tempmatchedisomerind(kk)];
                        end
                    end
                end
            end
            % Merge - If EThcD_cofrag is needed
            tempbestisomers_string = '';
            tempisomers_noncommonfound_numerical = zeros(1,length(isomers));
            tempisomers_numnoncommonfragfound_quant = zeros(1,length(isomers));
            if any(displaydataind_colCID)
                tempisomers_noncommonfound_numerical = tempisomers_noncommonfound_numerical + ...
                    CID_uniisomers_numerical + CID_sharedisomers_numerical;
                tempisomers_numnoncommonfragfound_quant = tempisomers_numnoncommonfragfound_quant + ...
                    CID_uniisomers_quant + CID_sharedisomers_quant;
            end
            if any(displaydataind_colHCD)
                tempisomers_noncommonfound_numerical = tempisomers_noncommonfound_numerical + ...
                    HCD_uniisomers_numerical + HCD_sharedisomers_numerical;
                tempisomers_numnoncommonfragfound_quant = tempisomers_numnoncommonfragfound_quant + ...
                    HCD_uniisomers_quant + HCD_sharedisomers_quant;
            end
            if any(displaydataind_colEThcD)
                tempisomers_noncommonfound_numerical = tempisomers_noncommonfound_numerical + ...
                    EThcD_uniisomers_numerical + EThcD_sharedisomers_numerical;
                tempisomers_numnoncommonfragfound_quant = tempisomers_numnoncommonfragfound_quant + ...
                    EThcD_uniisomers_quant + EThcD_sharedisomers_quant;
            end
            if any(displaydataind_colETciD)
                tempisomers_noncommonfound_numerical = tempisomers_noncommonfound_numerical + ...
                    ETciD_uniisomers_numerical + ETciD_sharedisomers_numerical;
                tempisomers_numnoncommonfragfound_quant = tempisomers_numnoncommonfragfound_quant + ...
                    ETciD_uniisomers_quant + ETciD_sharedisomers_quant;
            end
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
                    [uniidx,~,uniidxind] = unique(idx);
                    uniidx_meanscore = zeros(size(uniidx));
                    for jj = 1:length(uniidx)
                        uniidx_meanscore(jj) = mean(foundisomers_numnoncommonfragfound_quant(uniidxind == jj));
                    end
                    [~,highscoreind] = max(uniidx_meanscore);
                    tgtidx = uniidx(highscoreind);
                    tempbestisomers = foundisomers(idx == tgtidx);
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
            currentind = currentind + length(colnames);
            newprotid = str2num(result(displaydataind(ii,1)).ProteinID); %protein ID 1 5 266, result goes from 1 to 516 scans, ii = groups of scans #
            %             newprotid1 = newprotid(1:2);
            %             newprotid2 = newprotid(3:end);
            %             for kk = 1:length(newprotid2)
            %                 currentptmid = newprotid2(kk);
            %                 newprotid2(kk) = find(ismember(new_ptmseq,scoreintdata.ptmisomers{currentptmid,1}));
            %             end
            %             newprotid = [newprotid1,newprotid2];
            if any(displaydataind_colHCD)
                tempresult_HCD = result(displaydataind(ii,displaydataind_colHCD));
                tempresult_HCD.ProteinID = newprotid;
            end
            if any(displaydataind_colCID)
                tempresult_CID = result(displaydataind(ii,displaydataind_colCID));
                tempresult_CID.ProteinID = newprotid;
            end
            if any(displaydataind_colEThcD)
                tempresult_EThcD = result(displaydataind(ii,displaydataind_colEThcD));
                if strcmpi(tempresult_EThcD.Fragmode,'ethcd')
                    tempresult_EThcD.Fragmode = 'EThcD';
                end
                tempresult_EThcD.ProteinID = newprotid;
            end
            if any(displaydataind_colETciD)
                tempresult_ETciD = result(displaydataind(ii,displaydataind_colETciD));
                tempresult_ETciD.ProteinID = newprotid;
            end
            if any(displaydataind_colETD)
                tempresult_ETD = result(displaydataind(ii,displaydataind_colETD));
                tempresult_ETD.ProteinID = newprotid;
            end
            if currentind > sectionsize
                isomericresult = [isomericresult;tempresult];
                tempresult = cell(sectionsize,30);
                currentind = 1;
            end
            for kk = 1:length(fldnms)
                for ll = 1:length(colnames)
                    switch upper(colnames{ll})
                        case 'HCD'
                            tempresult{currentind + ll - 1,kk} = tempresult_HCD.(fldnms{kk});
                        case 'CID'
                            tempresult{currentind + ll - 1,kk} = tempresult_CID.(fldnms{kk});
                        case 'ETHCD'
                            tempresult{currentind + ll - 1,kk} = tempresult_EThcD.(fldnms{kk});
                        case 'ETD'
                            tempresult{currentind + ll - 1,kk} = tempresult_ETD.(fldnms{kk});
                        case 'ETCID'
                            tempresult{currentind + ll - 1,kk} = tempresult_ETciD.(fldnms{kk});
                    end
                end
            end
        else
            theofragstoptmseq = cell(1,length(theofrag_HCDsto));
            if any(displaydataind_colHCD)
                for jj = 1:length(theofrag_HCDsto)
                    theofragstoptmseq{jj} = theofrag_HCDsto{jj}(1).original;
                end
            elseif any(displaydataind_colCID)
                for jj = 1:length(theofrag_CIDsto)
                    theofragstoptmseq{jj} = theofrag_CIDsto{jj}(1).original;
                end
            elseif any(displaydataind_colEThcD)
                for jj = 1:length(theofrag_EThcD_cofragsto)
                    theofragstoptmseq{jj} = theofrag_EThcD_cofragsto{jj}(1).original;
                end
            elseif any(displaydataind_colETD)
                for jj = 1:length(theofrag_ETD_cofragsto)
                    theofragstoptmseq{jj} = theofrag_ETD_cofragsto{jj}(1).original;
                end
            elseif any(displaydataind_colETciD)
                for jj = 1:length(theofrag_ETciD_cofragsto)
                    theofragstoptmseq{jj} = theofrag_ETciD_cofragsto{jj}(1).original;
                end
            end
            tempbestisomers = strsplit(bestisomers,', ');
            newisomericglycan = setdiff(tempbestisomers,g_ori.struct);
            if any(displaydataind_colHCD)
                tempresult_HCD = result(displaydataind(ii,displaydataind_colHCD));
                tempresult_HCD.ProteinID = isomerProteinIDsto{1};
            end
            if any(displaydataind_colCID)
                tempresult_CID = result(displaydataind(ii,displaydataind_colCID));
                tempresult_CID.ProteinID = isomerProteinIDsto{1};
            end
            if any(displaydataind_colEThcD)
                tempresult_EThcD = result(displaydataind(ii,displaydataind_colEThcD));
                tempresult_EThcD.ProteinID = isomerProteinIDsto{1};
            end
            if any(displaydataind_colETD)
                tempresult_ETD = result(displaydataind(ii,displaydataind_colETD));
                tempresult_ETD.ProteinID = isomerProteinIDsto{1};
            end
            if any(displaydataind_colETciD)
                tempresult_ETciD = result(displaydataind(ii,displaydataind_colETciD));
                tempresult_ETciD.ProteinID = isomerProteinIDsto{1};
            end
            if ismember(g_ori.struct,tempbestisomers)
                currentind = currentind + length(colnames);
                if currentind > sectionsize
                    isomericresult = [isomericresult;tempresult];
                    tempresult = cell(sectionsize,30);
                    currentind = 1;
                end
                for kk = 1:length(fldnms)
                    for ll = 1:length(colnames)
                        switch upper(colnames{ll})
                            case 'HCD'
                                tempresult{currentind + ll - 1,kk} = tempresult_HCD.(fldnms{kk});
                            case 'CID'
                                tempresult{currentind + ll - 1,kk} = tempresult_CID.(fldnms{kk});
                            case 'ETHCD'
                                tempresult{currentind + ll - 1,kk} = tempresult_EThcD.(fldnms{kk});
                            case 'ETD'
                                tempresult{currentind + ll - 1,kk} = tempresult_ETD.(fldnms{kk});
                            case 'ETCID'
                                tempresult{currentind + ll - 1,kk} = tempresult_ETciD.(fldnms{kk});
                        end
                    end
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
                    %% HCD
                    if any(displaydataind_colHCD)
                        methodind = ismember(upper(scoreoptions.analyzefragmode),'HCD');
                        theofrag_HCD = theofrag_HCDsto{sgpind};
                        [spectrascores_HCD,~] = match1by1(thisnewglypepseq,tempresult_HCD.Theo,...
                            theofrag_HCD,isomericglypepprotid,scanHCD,...
                            tempresult_HCD.Expt,tempresult_HCD.Mono,tempresult_HCD.Retime,spectrumHCD,chgHCD,...
                            scoreoptions.ms2tol(methodind),scoreoptions.ms2tolunit{methodind},'HCD',scoreintdata.sliminput.fragnum,tempresult_HCD.MDiff,...
                            tempresult_HCD.Quant,tempresult_HCD.Protein,struct('maxlag',scoreoptions.maxlag,'selectpeak',[],'proceed_override',true));   % get the score for 1 spectrum vs 1 candidate for 1 fragmentation mode
                    end
                    %% CID
                    if any(displaydataind_colCID)
                        methodind = ismember(upper(scoreoptions.analyzefragmode),'CID');
                        theofrag_CID = theofrag_CIDsto{sgpind};
                        [spectrascores_CID,~] = match1by1(thisnewglypepseq,tempresult_CID.Theo,...
                            theofrag_CID,isomericglypepprotid,scanCID,...
                            tempresult_CID.Expt,tempresult_CID.Mono,tempresult_CID.Retime,spectrumCID,chgCID,...
                            scoreoptions.ms2tol(methodind),scoreoptions.ms2tolunit{methodind},'CID',scoreintdata.sliminput.fragnum,tempresult_CID.MDiff,...
                            tempresult_CID.Quant,tempresult_CID.Protein,struct('maxlag',scoreoptions.maxlag,'selectpeak',[],'proceed_override',true));
                    end
                    %% EThcD
                    if any(displaydataind_colEThcD)
                        methodind = ismember(upper(scoreoptions.analyzefragmode),'ETHCD');
                        theofrag_EThcD = theofrag_EThcD_cofragsto{sgpind};
                        [spectrascores_EThcD,~] = match1by1(thisnewglypepseq,tempresult_EThcD.Theo,...
                            theofrag_EThcD,isomericglypepprotid,scanEThcD,...
                            tempresult_EThcD.Expt,tempresult_EThcD.Mono,tempresult_EThcD.Retime,spectrumEThcD,chgEThcD,...
                            scoreoptions.ms2tol(methodind),scoreoptions.ms2tolunit{methodind},'EThcD',scoreintdata.sliminput.fragnum(methodind,:),tempresult_EThcD.MDiff,...
                            tempresult_EThcD.Quant,tempresult_EThcD.Protein,struct('maxlag',scoreoptions.maxlag,'selectpeak',[],'proceed_override',true));
                    end
                    %% ETD
                    if any(displaydataind_colETD)
                        methodind = ismember(upper(scoreoptions.analyzefragmode),'ETD');
                        theofrag_ETD = theofrag_ETD_cofragsto{sgpind};
                        [spectrascores_ETD,~] = match1by1(thisnewglypepseq,tempresult_ETD.Theo,...
                            theofrag_ETD,isomericglypepprotid,scanETD,...
                            tempresult_ETD.Expt,tempresult_ETD.Mono,tempresult_ETD.Retime,spectrumETD,chgETD,...
                            scoreoptions.ms2tol(methodind),scoreoptions.ms2tolunit{methodind},'ETD',scoreintdata.sliminput.fragnum(methodind,:),tempresult_ETD.MDiff,...
                            tempresult_ETD.Quant,tempresult_ETD.Protein,struct('maxlag',scoreoptions.maxlag,'selectpeak',[],'proceed_override',true));
                    end
                    %% ETciD
                    if any(displaydataind_colETciD)
                        methodind = ismember(upper(scoreoptions.analyzefragmode),'ETHCD');
                        theofrag_ETciD = theofrag_ETciD_cofragsto{sgpind};
                        [spectrascores_ETciD,~] = match1by1(thisnewglypepseq,tempresult_ETciD.Theo,...
                            theofrag_ETciD,isomericglypepprotid,scanETciD,...
                            tempresult_ETciD.Expt,tempresult_ETciD.Mono,tempresult_ETciD.Retime,spectrumETciD,chgETciD,...
                            scoreoptions.ms2tol(methodind),scoreoptions.ms2tolunit{methodind},'ETciD',scoreintdata.sliminput.fragnum(methodind,:),tempresult_ETciD.MDiff,...
                            tempresult_ETciD.Quant,tempresult_ETciD.Protein,struct('maxlag',scoreoptions.maxlag,'selectpeak',[],'proceed_override',true));
                    end
                    currentind = currentind + length(colnames);
                    if currentind > sectionsize
                        isomericresult = [isomericresult;tempresult];
                        tempresult = cell(sectionsize,30);
                        currentind = 1;
                    end
                    for kk = 1:length(fldnms)
                        for ll = 1:length(colnames)
                            switch upper(colnames{ll})
                                case 'HCD'
                                    tempresult{currentind + ll - 1,kk} = spectrascores_HCD.(fldnms{kk});
                                case 'CID'
                                    tempresult{currentind + ll - 1,kk} = spectrascores_CID.(fldnms{kk});
                                case 'ETHCD'
                                    tempresult{currentind + ll - 1,kk} = spectrascores_EThcD.(fldnms{kk});
                                case 'ETD'
                                    tempresult{currentind + ll - 1,kk} = spectrascores_ETD.(fldnms{kk});
                                case 'ETCID'
                                    tempresult{currentind + ll - 1,kk} = spectrascores_ETciD.(fldnms{kk});
                            end
                        end
                    end
                end
            end
        end
    end
end
isomericresult = [isomericresult;tempresult];
isomericresult = isomericresult(~cellfun(@isempty,isomericresult(:,1)),:);
ProteinIDcol = ismember(fldnms,'ProteinID');
for ii = 1:size(isomericresult,1)
    tempProteinID = isomericresult{ii,ProteinIDcol};
    if isnumeric(tempProteinID)
        isomericresult{ii,ProteinIDcol} = num2str(tempProteinID);
    end
end
isomericglypepresult = cell2struct(isomericresult,fldnms,2)';
%delete(wtbar);
end