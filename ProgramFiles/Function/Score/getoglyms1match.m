function [glypepcomp,glypepmass,glycancombi] = getoglyms1match(...
    pepseq,pepmass,precmass,PTM,glycancombi,glycancombimass,scoreoptions)
% pepseq: n x 1 cell array
% pepmass: n x 1 numerical array
% precmass: 1 number
% glycancombi: 1 x n cell array of numerical
glypepcomp = {};
glypepmass = [];
PTMuniind = PTM.ptmuniind;
PTMisvar = PTM.ptmisvar;
PTMisglycan = PTM.ptmisglycan;
PTMrealmass = PTM.ptmmass + PTM.ptmmassfix;
for ii = 1:length(pepseq)
    if ii == 28
        A=1;
    end
    thispepmass = pepmass(ii);
    thispepseq = pepseq{ii};
    pepseq_PTMseq = regexp(thispepseq,'\{\d+\}','match');  % find and locate PTMs
    ptmonpep_isglycan = false(size(pepseq_PTMseq));
    ptmonpep_isvar = false(size(pepseq_PTMseq));
    
    ptmuniind_onpep = zeros(size(pepseq_PTMseq));
    for jj = 1:length(pepseq_PTMseq)
        thisptmuniind_onpep = str2double(pepseq_PTMseq{jj}(2:end-1));  % 1, 2, 3,... represents PTM marker
        ptmuniind_onpep(jj) = thisptmuniind_onpep;  % looks like [2 2 1 1 2 3 ...]
        thisptmuniind_onpepPTMind = find(PTMuniind == thisptmuniind_onpep);
        ptmonpep_isglycan(jj) = PTMisglycan(thisptmuniind_onpepPTMind(1));
        ptmonpep_isvar(jj) = PTMisvar(thisptmuniind_onpepPTMind(1));
    end
    [uniptmuniind_onpep,~,uniptmuniind_onpepind] = unique(ptmuniind_onpep);
    ptmonpepstat = zeros(length(uniptmuniind_onpep),3);  % [uniind,#onpep,#avail]
    ptmonpepstat(:,1) = uniptmuniind_onpep;
    for jj = 1:max(uniptmuniind_onpepind)
        tempind = find(uniptmuniind_onpepind == jj);
        ptmonpepstat(jj,2) = length(tempind);
        ptmonpepstat(jj,3) = sum(PTMuniind == uniptmuniind_onpep(jj));
    end
    PTMmasstgt = precmass - thispepmass;
    fixedPTMcomp = [];
    fixedPTMmass = 0;
    ptmonpepstat_isvarind = true(size(ptmonpepstat,1),1);
    if any(~ptmonpep_isvar)  % fixed PTM on peptide
        fixedPTMonpepuniind = ptmuniind_onpep(~ptmonpep_isvar);
        for jj = 1:length(fixedPTMonpepuniind)
            fixedPTMpos = find(PTMuniind == fixedPTMonpepuniind(jj));
            PTMmasstgt = PTMmasstgt - PTMrealmass(fixedPTMpos);
            fixedPTMcomp = [fixedPTMcomp,fixedPTMpos];
        end
        [~,isvarind] = ismember(ptmonpepstat(:,1),fixedPTMonpepuniind);
        ptmonpepstat_isvarind(~isvarind) = false;
        fixedPTMmass = sum(PTMrealmass(fixedPTMcomp));
    end
    ptmonpepstat_isvarnonglyind = ptmuniind_onpep(~ptmonpep_isglycan & ptmonpep_isvar);
    ptmonpepstat_isvarnongly = ptmonpepstat(ismember(ptmonpepstat(:,1),ptmonpepstat_isvarnonglyind),:);
    ptmonpepstat_isvarglyind = ptmuniind_onpep(ptmonpep_isglycan & ptmonpep_isvar);
    ptmonpepstat_isvargly = ptmonpepstat(ismember(ptmonpepstat(:,1),ptmonpepstat_isvarglyind),:);
    % Above 2 looks like:
    % uniind    #onpep    #avail
    % [1                 2                1    
    %  2                  4              27]
    % varnongly: variable non glycan ptm, only 2 states: on/off
    if size(ptmonpepstat_isvargly,1) > 1
        error('Multiple glycosylation patterns detected, not supported in current version');
    end
    tempptmlayout_varnongly = cell(size(ptmonpepstat_isvarnongly,1),1);
    for jj = 1:size(ptmonpepstat_isvarnongly,1)
        tempptmlayout_varnongly{jj} = nmultichoosek(...
            find(PTMuniind == ptmonpepstat_isvarnongly(jj,1)),ptmonpepstat_isvarnongly(jj,2));
    end
    ptmlayout_varnongly = combigroups(tempptmlayout_varnongly);
    if isempty(ptmlayout_varnongly)
        ptmlayout_varnongly = 0;
        %ptmonpepstat_isvargly = [0,0,0];
    end
    for jj = 1:size(ptmlayout_varnongly,1)
        hitglycancombi = {};
        hitglycancombimass = [];
        numvarglysites = min(ptmonpepstat_isvargly(1,2),scoreoptions.maxOglyonpep);
        tempptmlayout_varnongly = ptmlayout_varnongly(jj,ptmlayout_varnongly(jj,:) > 0);
        tempptmlayout_varnonglymass = sum(PTMrealmass(tempptmlayout_varnongly));
        PTMmasstgt_wovarnongly = PTMmasstgt - tempptmlayout_varnonglymass;
        if numvarglysites <= length(glycancombi)
            if numvarglysites == 0 %if isempty(glycancombi{numvarglysites})
                tempglycancombi = getglycancombi(PTM,numvarglysites);
                glycancombi = tempglycancombi;
                tempglycancombirealmass = zeros(size(tempglycancombi,1),1);
                for kk = 1:size(tempglycancombi,1)
                    temptempglycancombi = 0; %tempglycancombi(kk,1)
                    temptempglycancombi = temptempglycancombi(temptempglycancombi > 0);
                    tempglycancombirealmass(kk) = sum(PTMrealmass(temptempglycancombi));
                end
                glycancombimass = tempglycancombirealmass; %{numvarglysites}
                tempglycancombi = glycancombi;
                tempglycancombimass = glycancombimass;
            else
                tempglycancombi = glycancombi{numvarglysites};
                tempglycancombimass = glycancombimass{numvarglysites};
            end
            switch upper(scoreoptions.ms1tolunit)
                case 'DA'
                    tgthitind = abs(tempglycancombimass - PTMmasstgt_wovarnongly) <= scoreoptions.ms1tol;
                case 'PPM'
                    tgthitind = abs(tempglycancombimass - PTMmasstgt_wovarnongly)/precmass * 1e6...
                        <= scoreoptions.ms1tol;
            end
            if any(tgthitind)
                tempglycancombihit = tempglycancombi(tgthitind,:);
                tempglycancombimasshit = tempglycancombimass(tgthitind);
                for kk = 1:size(tempglycancombihit,1)
                    temphitglycancombi = [ii,fixedPTMcomp,tempptmlayout_varnongly,...
                        tempglycancombihit(kk,:)];
                    temphitglycancombi = temphitglycancombi(temphitglycancombi > 0);
                    hitglycancombi = [hitglycancombi;temphitglycancombi];
                    hitglycancombimass = [hitglycancombimass;...
                        thispepmass + fixedPTMmass + tempptmlayout_varnonglymass + tempglycancombimasshit(kk)];
                end
            end
        else
            num1sthalf = ceil(numvarglysites/2);
            if (num1sthalf > length(glycancombi)) || isempty(glycancombi{num1sthalf})
                tempglycancombi = getglycancombi(PTM,num1sthalf);
                glycancombi{num1sthalf} = tempglycancombi;
                tempglycancombirealmass = zeros(size(tempglycancombi,1),1);
                for kk = 1:size(tempglycancombi,1)
                    temptempglycancombi = tempglycancombi(kk,1);
                    temptempglycancombi = temptempglycancombi(temptempglycancombi > 0);
                    tempglycancombirealmass(kk) = sum(PTMrealmass(temptempglycancombi));
                end
                glycancombimass{num1sthalf} = tempglycancombirealmass;
            end
            combi1sthalf = glycancombi{num1sthalf};
            mass1sthalf = glycancombimass{num1sthalf};
            num2ndhalf = numvarglysites - num1sthalf;
            if isempty(glycancombi{num2ndhalf})
                tempglycancombi = getglycancombi(PTM,num2ndhalf);
                glycancombi{num2ndhalf} = tempglycancombi;
                tempglycancombirealmass = zeros(size(tempglycancombi,1),1);
                for kk = 1:size(tempglycancombi,1)
                    temptempglycancombi = tempglycancombi(kk,1);
                    temptempglycancombi = temptempglycancombi(temptempglycancombi > 0);
                    tempglycancombirealmass(kk) = sum(PTMrealmass(temptempglycancombi));
                end
                glycancombimass{num2ndhalf} = tempglycancombirealmass;
            end
            combi2ndhalf = glycancombi{num2ndhalf};
            mass2ndhalf = glycancombimass{num2ndhalf};
            tgtmass2ndhalf = PTMmasstgt_wovarnongly - mass1sthalf;
            tgtmass2ndhalf_nonnegind = tgtmass2ndhalf >= 0;
            tgtmass2ndhalf_nonneg = tgtmass2ndhalf(tgtmass2ndhalf_nonnegind);
            combi1sthalf_nonneg = combi1sthalf(tgtmass2ndhalf_nonnegind);
            mass1sthalf_nonneg = mass1sthalf(tgtmass2ndhalf_nonnegind);
            
            for kk = 1:length(tgtmass2ndhalf_nonneg)
                switch upper(scoreoptions.ms1tolunit)
                    case 'DA'
                        tgthitind2ndhalf = abs(tgtmass2ndhalf_nonneg(kk) - mass2ndhalf) <= scoreoptions.ms1tol;
                    case 'PPM'
                        tgthitind2ndhalf = abs(tgtmass2ndhalf_nonneg(kk) - mass2ndhalf)/precmass * 1e6...
                            <= scoreoptions.ms1tol;
                end
                if any(tgthitind2ndhalf)
                    combi2ndhlafhit = combi2ndhalf(tgthitind2ndhalf,:);
                    temphitglycancombi = cell(size(combi2ndhlafhit));
                    temphitglycancombimass = zeros(size(combi2ndhlafhit));
                    for ll = 1:length(combi2ndhlafhit)
                        temptemphitglycancombi = [ii,fixedPTMcomp,tempptmlayout_varnongly,...
                            combi1sthalf_nonneg(kk,:),combi2ndhlafhit(ll,:)];
                        temphitglycancombi{ll} = temptemphitglycancombi(temptemphitglycancombi > 0);
                        temphitglycancombimass(ll) = thispepmass + fixedPTMmass + ...
                            mass1sthalf_nonneg(kk) + temphitglycancombimass(ll);
                    end
                    hitglycancombi = [hitglycancombi;temphitglycancombi];
                    hitglycancombimass = [hitglycancombimass;temphitglycancombimass];
                end
            end
        end
        glypepcomp = [glypepcomp;hitglycancombi];
        glypepmass = [glypepmass;hitglycancombimass];
    end
end
end