function result = analyzegp(matched_gpcomp,matched_gpmass,unisearchscannum,chargesto,spectrasto,...
    retimesto,exptmasssto,monomasssto,quantsto,allunisearchscannumind,...
    unisearchscannumindind,protbatch,FASTAbatch,protind,PTM,...
    agpoptions)
% ANALYZEGP: fragment glycopeptides, then match against experiment spectra
%
% Syntax:
% result = analyzegp(matched_gpcomp,matched_gpmass,unisearchscannum,chargesto,spectrasto,...
%     retimesto,exptmasssto,monomasssto,quantsto,allunisearchscannumind,...
%     unisearchscannumindind,protbatch,FASTAbatch,protind,PTM,...
%     agpoptions)
%
% Input:
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%     (Below the value of all "n"s and all "k"s are equal)
%     Variable name                  Format                                 Description
%     matched_gpcomp           n x 1 cell array of               Composition of matched candidates.
%                                                      strings.
%     matched_gpmass            n x 1 double                        Theoretical mass of matched candidates.
%     unisearchscannum          k x 1 double                         Scan numbers to be searced. These are unique values so
%                                                                                                     each scan number appears only once.
%     chargesto                           k x 1 double                        Charge state of the precursor ion.
%     spectrasto                          k x 1 cell array of               MS2 spectra.
%                                                      m x 2 double.
%     retimesto                            k x 1 double                        Retention time of precursor ion.
%     exptmasssto                       k x 1 double                        "Isolation window target" converted to massfrom m/z.
%                                                                                                      Here this value is treated as precursor ion mass.
%     monomasssto                    k x 1 double                        Monoisotopic mass of precursor ion.
%     quantsto                             k x 1 double                        AUC of precursor ion.
%     allunisearchscannumind n x m double                       Shows which unique scan number matches
%                                                                                                      "matched_gpcomp". Each row corresponds to one
%                                                                                                      gpcomp. For triggered experiment, each column refers
%                                                                                                      to one fragmentation method.
%     unisearchscannumindind 1 x n double                       Index of "allunisearchscannumind". This indicates which
%                                                                                                      part of the whole job this worker received.
%     protbatch                            m x 1 cell array of             Candidate peptide backbones.
%                                                      p x 1 array of strings.
%     FASTAbatch                        m x 1 cell array of              Protein FASTA names.
%                                                      strings.
%     protind                                 1 x m double                       Protein index.
%     PTM                                      Structure                             PTM information.
%     (See THEOPTMFRAG for details.)
%         ~.ptmseq                         1 x m cell array of              PTM sequence.
%                                                        strings
%         ~.ptmmass                      1 x m double                       PTM mass.
%         ~.ptmmassfix                  1 x m double                       Mass fix values when combining PTM with peptide
%                                                                                                       backbone.
%
%         ~.ptmtype                        1 x m double                       PTM type.
%         ~.ptmfragsto                   k x m cell array of              PTM fragments.
%                                                        1 x p structures
%         ~.ptmfragstublen           k x m cell array of              "Length" of each fragment.
%                                                        1 x p double
%         ~.ptmuniind                     m x 1 double                       PTM locator.
%     agpoptions
%         ~.ms1tolunit                    String                                   Tolerence of precursor ion matching (unit).
%         ~.ms2tol                            n x 1 double                       Tolerence of MS2 matching (value).
%         ~.ms2tolunit                     n x 1 cell array of              Tolerence of MS2 matching (value).
%                                                        strings
%         ~.fragnum                         n x m double                      Number of cleavages allowed for each fragmentation
%                                                                                                       method.
%         ~.fragmode                       n x 1 cell array of              Fragmentation modes.
%                                                        strings
%                                                  (Below: see SCOREALLSPECTRA for details.)
%         ~.selectpeak
%         ~.maxlag
%         ~.theofragopt_
%             monosacislabile
%         ~.theofragopt_
%             simultaneousfrag
%         ~.theofragopt_
%             addoxoniumion
%         ~.theofragopt_
%             maxstublen
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%
% Output:
% result: 1 x n structure. Analysis results. See MATCH1BY1 for detail.
% Result contains all parameters used for single spectrum scoring
% including: Scan [scannum]; Mono [monomass]; Theo [sgpmass]; Expt
% [exptmass]; Charge [charge]; SGP [sgp]; Fragmode = fragmethod; PeakLag;
% HtCenter; HtAvg; PercentIonMatch; Pvalue; DecoyRatio; Top10; SelectPeak;
% NpFrag [fragnum(1)]; NgFrag [fragnum(2)]; NmFrag [fragnum(3)]; Enscore; 
% DecoyEnscore; Retime; Quant [quant]; Glycov; Pepcov; Y0Y1Y2 = num2str([0,0,0]);
% ProteinID [proteinID]; Fracpepfragmatch; Fracglyfragmatch; MDiff [mdiff]; Protein [proteinname];
% 
%
% Note:
% N/A
%
% Example:
% N/A
%
% Children function:
% N/A
% See also:
% 

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 04/01/2020

blocksize = 100000;
fldnms = ResultItem.allresultitems;

unisearchscannumind = allunisearchscannumind(unisearchscannumindind,:);
glypepcomplength = cellfun(@length,matched_gpcomp);
glypepcomp = zeros(length(matched_gpcomp),max(glypepcomplength));
for ii = 1:length(matched_gpcomp)
    glypepcomp(ii,1:glypepcomplength(ii)) = matched_gpcomp{ii};
end
result = [];
[unipepcomp,~,unipepind] = unique(glypepcomp(:,1:2),'stable','rows');
[uniglypepcomp,~,uniglypepind] = unique(glypepcomp,'stable','rows');
[unipepfragnum,~,unipepfragnumind] = unique(agpoptions.fragnum(:,1),'stable');
pepfragsto = cell(size(unipepcomp,1),length(unipepfragnum));
%% CAUTION! THE "FRAGPEP" IS USING EMPTY STRING AS FRAGMODE INPUT
% THIS IS ONLY A TEMPORARY SOLUTION
for ii = 1:length(unipepfragnum)
    for jj = 1:size(unipepcomp,1)
        pepseq = protbatch{protind == unipepcomp(jj,1)}{unipepcomp(jj,2)};
        pepfragsto{jj,ii} = fragpep(pepseq,unipepfragnum(ii),'bcyzi','');
    end
end

% CID with #monosac < 4 requires npFrag be changed to >= 1
if ismember(1,unipepfragnum)
    specialpepfragsto = pepfragsto(:,ismember(unipepfragnum,1));
else
    specialpepfragsto = cell(length(pepfragsto),1);
    for ii = 1:length(pepfragsto)
        pepseq = protbatch{protind == unipepcomp(ii,1)}{unipepcomp(ii,2)};
        specialpepfragsto{ii,1} = fragpep(pepseq,1,'bcyzi','');
    end
end

searchfragmodes = agpoptions.fragmode;
fragopt.stublen = agpoptions.theofragopt_maxstublen;
fragopt.mode = 2;
fragopt.fragmode = '';
fragopt.monosacislabile = false;
fragopt.simultaneousfrag = false;
fragopt.addoxoniumion = false;

currentind = blocksize + 1;
tempaccuresult = [];
stopsignal = false(size(matched_gpcomp,1),1);
for ii = 1:length(searchfragmodes)
    pepfragcolumnnum = unipepfragnumind(ii);
    fragopt.fragmode = searchfragmodes{ii};
    fragopt.monosacislabile = agpoptions.theofragopt_monosacislabile(ii);
    fragopt.simultaneousfrag = agpoptions.theofragopt_simultaneousfrag(ii);
    fragopt.addoxoniumion = agpoptions.theofragopt_addoxoniumion(ii);
    ms2tol = agpoptions.ms2tol(ii);
    ms2tolunit = agpoptions.ms2tolunit{ii};
    thisfragmode = searchfragmodes{ii};
    for jj = 1:max(uniglypepind)
        thisgpcomp_uniglypepind = find(uniglypepind == jj);
        if any(~stopsignal(thisgpcomp_uniglypepind))
            thisfragnum = agpoptions.fragnum(ii,:);
            thisuniglypepcomp = uniglypepcomp(jj,:);
            thisuniglypepcomp(thisuniglypepcomp == 0) = [];
            thisunicompind = find(uniglypepind == jj,1);
            pepfragrownum = unipepind(thisunicompind);
            if length(thisuniglypepcomp) == 2
                temptheofrag = pepfragsto{pepfragrownum,pepfragcolumnnum};
                if strcmpi(searchfragmodes{ii},'CID')  % CID special: too few monosac in CID will
                    temptheofrag = specialpepfragsto{pepfragrownum,1};
                    thisfragnum(1) = max(1,thisfragnum(1));
                end
            else
                temppepfrag = pepfragsto{pepfragrownum,pepfragcolumnnum};
                if strcmpi(searchfragmodes{ii},'CID')  % CID special: too few monosac in CID will
                    % cause peptide fragmentation
                    nummonosac = sum(PTM.ptmsize(thisuniglypepcomp(3:end)));
                    if nummonosac < 4
                        temppepfrag = specialpepfragsto{pepfragrownum,1};
                        thisfragnum(1) = max(1,thisfragnum(1));
                    end
                end
                temptheofrag = combitheofrag(temppepfrag,...
                    PTM.ptmfragsto(ii,:),PTM.ptmfragstublen(ii,:),...
                    protbatch{protind == thisuniglypepcomp(1)}{thisuniglypepcomp(2)},PTM.ptmseq,...
                    PTM.ptmtype,PTM.ptmmass,thisuniglypepcomp,thisfragnum,fragopt);
            end
            tempunisearchscannumind = unisearchscannumind(uniglypepind == jj,ii);
            tempunisearchscannumind = tempunisearchscannumind(tempunisearchscannumind > 0);
            if any(tempunisearchscannumind)
                tempspectra = spectrasto(tempunisearchscannumind);
                tempcharge = chargesto(tempunisearchscannumind);
                tempscannum = unisearchscannum(tempunisearchscannumind);
                tempexptmass = exptmasssto(tempunisearchscannumind);
                tempmonomass = monomasssto(tempunisearchscannumind);
                tempretime = retimesto(tempunisearchscannumind);
                tempquant = quantsto(tempunisearchscannumind);
                sgpmass = matched_gpmass(find(uniglypepind == jj,1));
                proteinname = FASTAbatch{protind == thisuniglypepcomp(1)};
                for kk = 1:length(tempspectra)
                    if ~stopsignal(thisgpcomp_uniglypepind(kk))
                        sgp = temptheofrag(1).original;
                        thistempmonomass = tempmonomass(kk);
                        switch upper(agpoptions.ms1tolunit)
                            case 'DA'
                                mdiff = thistempmonomass - sgpmass;
                            case 'PPM'
                                mdiff = (thistempmonomass - sgpmass)/thistempmonomass*1e6;
                        end
                        [tempresult,tempstopsignal] = match1by1(sgp,sgpmass,temptheofrag,num2str(thisuniglypepcomp),tempscannum(kk),...
                            tempexptmass(kk),thistempmonomass,tempretime(kk),tempspectra{kk},tempcharge(kk),...
                            ms2tol,ms2tolunit,thisfragmode,thisfragnum,mdiff,...
                            tempquant(kk),proteinname,agpoptions);
                        if tempstopsignal
                            stopsignal(thisgpcomp_uniglypepind(kk)) = true;
                        end
                        currentind = currentind - 1;
                        if currentind == 0
                            result = [result,tempaccuresult];
                            currentind = blocksize;
                            tempaccuresult = [];
                        end
                        for mm = 1:length(fldnms)
                            tempaccuresult(currentind).(fldnms{mm}) = tempresult.(fldnms{mm});
                        end
                    end
                end
            end
        end
    end
end
tempaccuresult = tempaccuresult(~cellfun(@isempty,{tempaccuresult.Scan}));
result = [result,tempaccuresult];
end