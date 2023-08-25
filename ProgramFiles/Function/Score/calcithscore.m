function result = calcithscore(SpectraA,theofragmz,charge,tol,tolunit,mode,options)
% CALCITHSCORE: match theoretical fragments to experimental data.
% Quantifies top10, XCorr parameter values, ion match
% 
% Syntax:
% result = calcithscore(SpectraA,theofragmz,charge,tol,tolunit,mode,...
%     options)
% 
% Input:
% SpectraA: n x 2 numerical array. (Denoised) Experiment spectrum.
% theofragmz: n x 1 numerical array. The m/z values of theoretical
% fragments. Duplicated values are OK.
% charge: integer. Charge of parent ion.
% tol: number. Tolerence value.
% tolunit: string. Tolerence unit, "Da" or "ppm".
% mode: integer. Function working mode/cases. 
%     Case 1: all results calculated. Used when we need full scoring
%     Case 2: ionmatchindex and peakmatchindex. For browserGUI for red/green plot display
%     Case 3: ionmatchindex data only. For decoy calculations
% options: structure. Must contain 2 fields:
%     maxlag: double. Maximum lag in cross-correlation calculation.
%     selectpeak: 1 x n double. Specific m/z values to search for.
% 
% Output:
% result: structure. Results from Matching. Fields are:
% peakLag, htCenter, htAvg: [peakLag,htCenter,htAvg] = CrossCorr(SpectraA,SpectraB,50,tolunit);
% found: number of matched theoretical fragments that's within spectrum m/z range.
% search: number of searched theoretical fragments that's within spectrum m/z range.
% percentIonMatch: percentage of matched vs. searched theoretical
% fragments (0~100).
% Top10: How many of the 10 tallest peaks in spectrum were matched.
% ionmatchindex: 1 x n cell array. Which theoretical fragment is matched &
% what are those peaks.
% peakmatchindex: 1 x n cell array. Which peak is matched & who are those
% theoretical fragments.
% In mode 2 & 3 only part of the results will be calculated. Mode 2 returns
% "ionmatchindex", "peakmatchindex" and "selectpeakismatched". Mode 3
% returns "ionmatchindex" only.
% 
% Note:
% In mode 1 & 2, if a peak is matched and if the tolerence unit is "ppm",
% based on it's proposed charge state, program will try to look for its
% isotopic brethrens (within tolerence). Those isotopic peaks found will be
% considered as "matched" as well.
%
% Example:
% load "testdata_calcithscore.mat" and run command:
% result = calcithscore(SpectraA,theofragmz,charge,tol,tolunit,mode,options)
%
% Children function: 
% N/A

%% PREPARATION
switch lower(tolunit)
    case 'da'
        mzmin = min(SpectraA(:,1)) - tol;
        mzmax = max(SpectraA(:,1)) + tol;
    case 'ppm'
        mzmin = min(SpectraA(:,1)) * (1 - tol*1e-6);
        mzmax = max(SpectraA(:,1)) * (1 + tol*1e-6);
end
alltheofragmz = zeros(length(theofragmz),charge);
alltheofragchg = ones(length(theofragmz),charge);
for ii = 1:charge
    alltheofragmz(:,ii) = (theofragmz(:) - 1.007825032)/ii + 1.007825032;
    % theofrag are of chg = 1, this step calculates all possible m/z's at
    % different charge
    alltheofragchg(:,ii) = alltheofragchg(:,ii) * ii;
end
theofragind = repmat((1:length(theofragmz))',1,charge);  % they are the same ion even though chg be different
maxlag = options.maxlag;
selectpeak = options.selectpeak;

theofragmzfilter = alltheofragmz >= mzmin & alltheofragmz <= mzmax;  % delete m/z too high/low
filteredtheofragmz = alltheofragmz(theofragmzfilter);
filteredtheofragmz = filteredtheofragmz(:);
[filteredtheofragmz,uniind,~] = unique(filteredtheofragmz);  % this is to reduce computaion load
filteredtheofragchg = alltheofragchg(theofragmzfilter);
filteredtheofragchg = filteredtheofragchg(uniind);
tempfilteredtheofragind = theofragind(theofragmzfilter);
filteredtheofragind = tempfilteredtheofragind(uniind);
[~,~,isomerfragind] = unique(theofragmz,'stable');
% since we use unique m/z for calculation, matched m/z may refer to
% different fragments. the above step is to recover the association
% filteredtheofragmz: list of theoretical fragments in candidate
% filteredtheofragchg: charge of the theoretical fragments in candidate
% filteredtheofragind: index values for the theoretical fragments 
%% MATCH
switch mode
    case 1  % regular matching + return matching results
        % case 2 and case 3 are simplified versions of the following
        % procedure
        
        % Since only fragmz comes in, returned are indices
        ionismatched = cell(size(theofragmz));
        % "ionismatched" refers to theo fragments
        peakismatched = cell(size(SpectraA,1),1);
        % "peakismatched" refers to expt spectrum
        peakistheomatch = zeros(size(SpectraA,1),1);
        % "peakistheomatch" is similar to "peakismatched", but don't
        % consider isotopes
        selectpeakismatched = false(size(selectpeak));
        SpectraB = zeros(length(theofragmz),2);
        switch lower(tolunit)
            case 'ppm'
                for ii = 1:length(filteredtheofragmz)
                    thisfmz = filteredtheofragmz(ii);
                    thisfchg = filteredtheofragchg(ii);
                    ismatched = abs(SpectraA(:,1) - thisfmz)/thisfmz * 1e6 <= tol;
                    if any(ismatched)
                        whoismatched = find(ismatched);
                        ionismatched{filteredtheofragind(ii)} = [ionismatched{filteredtheofragind(ii)};whoismatched];
                        peakistheomatch(ismatched) = 1;
                        % isotope scenario
                        for jj = 1:length(whoismatched)
                            peakismatched{whoismatched(jj)} = [peakismatched{whoismatched(jj)};...
                                [filteredtheofragind(ii),filteredtheofragchg(ii)]];
                            keepisoMatch = true;  % the search stops as soon as proposed isotope is not found
                            isoplus = 0;
                            while keepisoMatch
                                isoplus = isoplus + 1;  % adding neutrons...
                                iso_thisfmz = thisfmz+1.007825032*isoplus/thisfchg;  % m/z after considering isotopes
                                iso_ismatched = abs(SpectraA(:,1) - iso_thisfmz)/iso_thisfmz * 1e6 <= tol;
                                if any(iso_ismatched)
                                    whoiso_ismatched = find(iso_ismatched);
                                    for k = 1:length(whoiso_ismatched)
                                        peakismatched{whoiso_ismatched(k)} = [peakismatched{whoiso_ismatched(k)};...
                                            [filteredtheofragind(ii),filteredtheofragchg(ii)]];
                                    end
                                else
                                    keepisoMatch = false;
                                end
                            end
                        end
                    end
                end
                for ii = 1:length(selectpeak)
                    thisfmz = selectpeak(ii);
                    ismatched = abs(SpectraA(:,1) - thisfmz)/thisfmz * 1e6 <= tol;
                    if any(ismatched)
                        selectpeakismatched(ii) = true;
                    end
                end
                for ii = 1:length(theofragmz)
                    if ~isempty(ionismatched{ii})
                        matchedpk = SpectraA(ionismatched{ii},:);                        
                        [~,ind] = max(matchedpk(:,2));
                        SpectraB(ii,1) = matchedpk(ind,1);
                    end
                end
            case 'da'
                for ii = 1:length(filteredtheofragmz)
                    thisfmz = filteredtheofragmz(ii);
                    ismatched = abs(SpectraA(:,1) - thisfmz) <= tol;
                    if any(ismatched)
                        whoismatched = find(ismatched);
                        ionismatched{filteredtheofragind(ii)} = [ionismatched{filteredtheofragind(ii)};whoismatched];
                        for jj = 1:length(whoismatched)
                            peakismatched{whoismatched(jj)} = [peakismatched{whoismatched(jj)};...
                                [filteredtheofragind(ii),filteredtheofragchg(ii)]];
                        end
                    end
                end
                for ii = 1:length(selectpeak)
                    thisfmz = selectpeak(ii);
                    ismatched = abs(SpectraA(:,1) - thisfmz) <= tol;
                    if any(ismatched)
                        selectpeakismatched(ii) = true;
                    end
                end
                for ii = 1:size(SpectraB,1)
                    if ~isempty(ionismatched{ii})
                        matchedpk = SpectraA(ionismatched{ii},:);                        
                        [~,ind] = max(matchedpk(:,2));
                        SpectraB(ii,1) = matchedpk(ind,1);
                    end
                end
        end
        for ii = 1:length(ionismatched)
            tempmatchedpks = ionismatched{ii};
            if ~isempty(tempmatchedpks)
                ionismatched(isomerfragind == isomerfragind(ii)) = {tempmatchedpks};
            end
        end
        found = sum(~cellfun(@isempty,ionismatched));
        search = length(ionismatched);
        percentionmatched = found / search * 100;
        Top10 = calculatetop10match(SpectraA,peakismatched,2,400);
        SpectraB = SpectraB(SpectraB(:,1)>0,:);
        SpectraB(:,2) = .1*max(SpectraA(:,2));
        if size(SpectraB,1)<3    % no need to do XCorr
            peakLag  = -50;
            htCenter = 0;
            htAvg    = 0;
        else
            [peakLag,htCenter,htAvg]=CrossCorr(SpectraA,SpectraB,maxlag,tolunit);
        end   
        result.peakLag = peakLag;  %a number (1st output) from CrossCorr function
        result.htCenter = htCenter; %a number (2nd output) from CrossCorr function
        result.htAvg = htAvg; %a number (3rd output) from CrossCorr function
        result.found = found; % the number of theoretical peaks that are matched
        result.search = search; % total number of theoretical peaks after filtering range
        result.percentIonMatch = percentionmatched;
        result.Top10 = Top10;  % a number from calculatetop10match stating how many of the Top10 peaks are matched
        result.ionmatchindex = ionismatched;  % see case 2 below for description
        result.peakmatchindex = peakismatched; % see case 2 below for description
        result.selectpeakismatched = selectpeakismatched; % see case 2 below for description
    case 2  % faster - return only matched index (peak & ion), for viewer
        ionismatched = cell(size(theofragmz));
        peakismatched = cell(size(SpectraA,1),1);
        selectpeakismatched = false(size(selectpeak));
        switch lower(tolunit)
            case 'ppm'
                for ii = 1:length(filteredtheofragmz)
                    thisfmz = filteredtheofragmz(ii);
                    thisfchg = filteredtheofragchg(ii);
                    ismatched = abs(SpectraA(:,1) - thisfmz)/thisfmz * 1e6 <= tol;
                    if any(ismatched)
                        whoismatched = find(ismatched);
                        ionismatched{filteredtheofragind(ii)} = [ionismatched{filteredtheofragind(ii)};whoismatched];
                        % isotope scenario
                        for jj = 1:length(whoismatched)
                            peakismatched{whoismatched(jj)} = [peakismatched{whoismatched(jj)};...
                                [filteredtheofragind(ii),filteredtheofragchg(ii)]];
                            keepisoMatch = true;
                            isoplus = 0;
                            while keepisoMatch
                                isoplus = isoplus + 1;
                                iso_thisfmz = thisfmz+1.007825032*isoplus/thisfchg;
                                iso_ismatched = abs(SpectraA(:,1) - iso_thisfmz)/iso_thisfmz * 1e6 <= tol;
                                if any(iso_ismatched)
                                    whoiso_ismatched = find(iso_ismatched);
                                    for k = 1:length(whoiso_ismatched)
                                        peakismatched{whoiso_ismatched(k)} = [peakismatched{whoiso_ismatched(k)};...
                                            [filteredtheofragind(ii),filteredtheofragchg(ii)]];
                                    end
                                    ionismatched{filteredtheofragind(ii)} = [ionismatched{filteredtheofragind(ii)};whoiso_ismatched];
                                else
                                    keepisoMatch = false;
                                end
                            end
                        end
                    end
                end
                for ii = 1:length(selectpeak)
                    thisfmz = selectpeak(ii);
                    ismatched = abs(SpectraA(:,1) - thisfmz)/thisfmz * 1e6 <= tol;
                    if any(ismatched)
                        selectpeakismatched(ii) = true;
                    end
                end
            case 'da'
                for ii = 1:length(filteredtheofragmz)
                    thisfmz = filteredtheofragmz(ii);
                    ismatched = abs(SpectraA(:,1) - thisfmz) <= tol;
                    if any(ismatched)
                        whoismatched = find(ismatched);
                        ionismatched{filteredtheofragind(ii)} = [ionismatched{filteredtheofragind(ii)};whoismatched];
                        for jj = 1:length(whoismatched)
                            peakismatched{whoismatched(jj)} = [peakismatched{whoismatched(jj)};...
                                [filteredtheofragind(ii),filteredtheofragchg(ii)]];
                        end
                    end
                end
                for ii = 1:length(selectpeak)
                    thisfmz = selectpeak(ii);
                    ismatched = abs(SpectraA(:,1) - thisfmz)/thisfmz * 1e6 <= tol;
                    if any(ismatched)
                        selectpeakismatched(ii) = true;
                    end
                end
        end
        for ii = 1:length(ionismatched)
            tempmatchedpks = ionismatched{ii};
            if ~isempty(tempmatchedpks)
                ionismatched(isomerfragind == isomerfragind(ii)) = {tempmatchedpks};
            end
        end
        result.ionmatchindex = ionismatched;
        result.peakmatchindex = peakismatched;
        result.selectpeakismatched = selectpeakismatched;
        % output is an array of 0s and 1s, with 1s indicating candidate fragments that are matched
        % theoretical array: [900, 915, 930, 930.1, 1000]
        % expt data: [899.9, 900.1 930.05, 1000]
        % result.ionmatchindex: cell array with matches:
        % {1,2},{},{3},{3},{4}-- 1 element for each theo fragment
        % result.peakmatchindex: cell array describing theo peak and charge state that matches each expt peak
        % {1,1+}{1,1+}{3,1+;4,1+}{5,1+} -- 1 element for each expt fragment 
        % selectpeak corresponds to the specific m/z values we are looking
        % for-- see row 21
    case 3  % fastest, only ionmatchindex, for decoy calculation use
        ionismatched = zeros(size(theofragmz));
        switch lower(tolunit)
            case 'ppm'
                for ii = 1:length(filteredtheofragmz)
                    thisfmz = filteredtheofragmz(ii);
                    ismatched = abs(SpectraA(:,1) - thisfmz)/thisfmz * 1e6 <= tol;
                    if any(ismatched)
                        ionismatched(filteredtheofragind(ii)) = 1;
                    end
                end
            case 'da'
                for ii = 1:length(filteredtheofragmz)
                    thisfmz = filteredtheofragmz(ii);
                    ismatched = abs(SpectraA(:,1) - thisfmz) <= tol;
                    if any(ismatched)
                        ionismatched(filteredtheofragind(ii)) = 1;
                    end
                end
        end
        for ii = 1:length(ionismatched)
            tempmatchedpks = ionismatched(ii);
            if tempmatchedpks > 0
                ionismatched(isomerfragind == isomerfragind(ii)) = 1;
            end
        end
        result.ionmatchindex = ionismatched;  
        % output is an array of 0s and 1s, with 1s indicating candidate fragments that are matched
        % for theoretical array [900, 915, 930, 930.1, 1000] and expt data
        % [899.9, 900.1 930.05, 1000], the output would be [1 0 1 1 1]
        % corresponding to the candidate that we are trying to match
end
end