function [defaultmimz,theomassdist,decision] = getdefaultmonoisowdist(precscannum,mz,chg,...
    allscannum,allcharge,allretime,allmslevel,allspectra,defaultmonoisodist,dist_range,...
    distyn,chnos)
%GETDEFAULTMONOISOWDIST: Performing Averagine model
%     calculations to correctify the precursor m/z values for a single scan
%     based on a supplied isotopic distribution
%
% Syntax:
% [defaultmimz,theomassdist] = getdefaultmonoisowdist(scannum,mz,chg,...
%     allscannum,allcharge,allretime,allmslevel,allspectra,defaultmonoisodist,dist_range,...
%     distyn,chnos)
%
% Input:
% scannum: double. The scan number to be analyzed.
% mz: double. The precursor ion's m/z (before correction)
% chg: double. Precursor ion's charge state.
% allscannum: n x 1 numerical array. All scan numbers.
% allcharge: n x 1 numerical array. All scan's charge state.
% allretime: n x 1 numerical array. All scan's retention time.
% allmslevel: n x 1 numerical array. All scan's MS level.
% allspectra: n x 1 cell array of m x 2 numerical array. All spectra
% defaultmonoisodist: n x 1 cell array of m x 2 numerical array. Isotope
%     distribution at different masses. The m x 2 numerical array is the
%     calculated isotope distribution whose monoisotopic mass is of a
%     specific mass, e.g. 4000 Da. In this array, 1st column is the mass,
%     2nd column is relative abundance.
% dist_range: double. The gap between adjacent distributions. In the
%     default model distr_range = 1Da.
% distyn: logical. If is true, only the isotopic peaks with abundance >=
%     10%
% chnos: 1 x 5 numerical array. User custom averagine model.
%
% Output:
%   defaultmimz: monoisotopic m/z that is consensus resilt
%   theomassdist: theoretical mass distributin
%   decision: rationale for the decision made ('A/B/C')
%
% Examples:
%   [defaultmimz,theomassdist] = getdefaultmonoisowdist(1569,759.1735,3,msdata.scannum,msdata.retime,msdata.mslevel,msdata.spectra,default_iso_dist,1,1,[])
%
%See also: PREPROCESSGUI, DEFAULTMONOISO, QUANTITATIVEANALYSIS
%
theomassdist = [];
defaultmimz = mz;
decision = 'N/A';
SearchConsecusiveMS1Scans = 6;
% If this value is set to ZERO, MS1 accumulation method is +/- 15 sec,
% Otherwise (value is integer and > 0), MS1 accumulation method is parent
% scan plus next N scans (SearchConsecusiveMS1Scans = N)

mass = mz * chg - 1.0078246*chg; % mass calculation
idx = floor(mass); % find corresponding theoretical distribution corresponding to mass
mslvl1ind = allmslevel == 1;
mslvl1ind2 = find(mslvl1ind);
if idx ~= 0 % if the precursor m/z is 0
    x = (precscannum == allscannum); % scan number index
    if SearchConsecusiveMS1Scans == 0
        tempretime = allretime(x); % retention time of scan number
        temprange = (abs(allretime-tempretime) <= .25) & mslvl1ind; % retention time window for all MS1 scans
    else
        thisprecscannumind = find(x);
        ms1scanstart = find(mslvl1ind2 == thisprecscannumind);
        temprange = mslvl1ind2(ms1scanstart:min(length(mslvl1ind2),ms1scanstart + SearchConsecusiveMS1Scans));
    end
    spectrain=allspectra(temprange); % MS1 spectra in retention time window
    spectrasize = cellfun(@size,spectrain,num2cell(ones(size(spectrain)))); % calculate size of each of these spectra
    spectra = zeros(sum(spectrasize),2); % create empty element with size equal to the sum of the elements in all spectra
    ind = 1; % dummy index
    for i = 1:length(spectrasize) % merge all input spectra in a single variable.
        spectra(ind:ind + spectrasize(i) - 1,:) = double(spectrain{i});
        ind = ind + spectrasize(i);
    end

    spectrumout = spectra(spectra(:,2) > 0,:); % keep only spectra that is not empty
    spectrumsectionind = ((spectrumout(:,1) - mz) <= 5/chg) & ((spectrumout(:,1) - mz) >= -3.1/chg); % section of spectra that is within tolerance of precursor m/z
    spectrumsection = spectrumout(spectrumsectionind,:);  % consider chg >= 2, 3 Da off monoiso max
    if idx <= defaultmonoisodist{end,1}(1,1) % if the precursor m/z is within the theoretical distribution
        idx_alt = round(idx/dist_range); % find index of mass based on distribution range
        if idx_alt == 0 % if rounding results in 0
            idx_alt = 1; % set value to 1
        elseif idx_alt > length(defaultmonoisodist) % if rounding results in value greater than length of theoretical distribution
            idx_alt = length(defaultmonoisodist); % set to max
        end
        theoisodist = defaultmonoisodist{idx_alt,1}; % select theoretical distribution
    else
        % [0.024321959339170,0.476559048284606,0.0014616980 46926,0.012941719677186,0] GLYCAN CHNOS VALUES
        % [0.0411425801242450,0.0667414095012345,0.00760993965077790,0.0203467465280654,0.000219515203030292] GLYCOPEPTIDE CHNOS VALUES
        if ~isempty(chnos)
            avggp = chnos;
        else
            avggp = [0.0411425801242450,0.0667414095012345,0.00760993965077790,0.0203467465280654,0.000219515203030292]; % [C H N O S]
        end
        [theoisodist,~,~] = isotopicdist(mass*avggp,'ShowPlot',0);
    end
    mzfix = (mass - floor(mass)) / chg;
    theoisodist(:,1) = theoisodist(:,1)/chg + 1.0078246 + mzfix;
    theoisostep = mean(diff(theoisodist(:,1)));
    expandtheoiso = flip(theoisodist(1,1)-theoisostep:-theoisostep:min(spectrumsection(:,1)));
    headlength = length(expandtheoiso);
    expandtheoiso = [expandtheoiso(:),zeros(length(expandtheoiso),1)];
    expandtheoiso = [expandtheoiso;theoisodist];
    etisotail = theoisodist(end,1)+theoisostep:theoisostep:max(spectrumsection(:,1))+1;
    expandtheoiso = [expandtheoiso;[etisotail(:),zeros(length(etisotail),1)]];
    expandtheoiso = unique(expandtheoiso,'rows');
    % Collecting all peaks in expandtheoiso about 10%
    if distyn == 1
        totIC = sum(expandtheoiso(:,2));
        theomassdist = [expandtheoiso(:,1),(expandtheoiso(:,2)/totIC)];
        intidx = any(theomassdist(:,2) < 0.1,2);
        theomassdist(intidx,:) = [];
    end
    expsmpl = zeros(size(expandtheoiso));
    expsmpl(:,1) = expandtheoiso(:,1);
    % expsmpl(:,1) = expsmpl(:,1) + mzfix;
    for j = 1:size(expsmpl,1)
        tempind = abs(spectrumsection(:,1) - expsmpl(j,1)) <= 0.0125;
        %         expsmpl(j,2) = sum(spectrumsection(tempind,2));
        if any(tempind)
            expsmpl(j,2) = max(spectrumsection(tempind,2));
        end
    end
    %     expsmplX = expsmpl;
    %     expsmplX(expsmplX(:,2) < 0.1 * max(expsmplX(:,2)),2) = 0;
    %     [~,ind] = sort(expsmplX(:,2),'descend');
    %     expsmplXmz = diff(sort(expsmplX(ind(1:4),1)));
    %     if max(expsmplXmz) - min(expsmplXmz) > 0.1
    %         go_ahead = false;
    %     else
    %         go_ahead = true;
    %     end
    go_ahead = true;
    if go_ahead
        % method A -  CLASSIC XCORR (FOR SIGNAL STRENGTH >= 10%)
        expsmplA = expsmpl;
        if mass<3000
            expsmplA(expsmpl(:,2) < 0.1 * max(expsmpl(:,2)),2) = 0;
        else
            expsmplA(expsmpl(:,2) < 0.05 * max(expsmpl(:,2)),2) = 0;
        end
        expandtheoisoA = expandtheoiso;
        expandtheoisoA(expandtheoiso(:,2) < 0.05 * max(expandtheoiso(:,2)),2) = 0;
        [a,b] = xcorr(expsmplA(:,2),expandtheoisoA(:,2),size(expandtheoisoA,1),'normalized');
        [~,ind] = max(a);
        lagA = b(ind);
        if headlength + 1 + lagA <= 0
            lagA = -headlength;
        end
        whichpeak = expsmplA(headlength + 1 + lagA,2);
        if whichpeak > 0
            defaultmimzA = theoisodist(1,1) + lagA * theoisostep;
        else  % fix by looking at higher m/z
            defaultmimzA = mz;
        end

        % method B - JUST TAKE THE LOWEST MZ
        expsmplB = expsmpl;
        if mass<3000
            expsmplB(expsmpl(:,2) < 0.1 * max(expsmpl(:,2)),:) = [];
        else
            expsmplB(expsmpl(:,2) < 0.03 * max(expsmpl(:,2)),:) = [];
        end
        defaultmimzB = expsmplB(1,1);

        % method C - TOP 4
        expsmplC = expsmpl;
        [~,ind] = sort(expsmplC(:,2),'descend');
        expsmplC(setdiff(1:size(expsmplC,1),ind(1:4)),2) = 0;
        [a,b] = xcorr(expsmplC(:,2),expandtheoiso(:,2),size(expandtheoiso,1),'normalized');
        [~,ind] = max(a);
        lagC = b(ind);
        if headlength + 1 + lagC <= 0
            lagC = -headlength;
        end
        whichpeak = expsmplA(headlength + 1 + lagC,2);
        if whichpeak > 0
            defaultmimzC = theoisodist(1,1) + lagC * theoisostep;
        else  % fix by looking at higher m/z
            defaultmimzC = mz;
        end
        agreeAB = abs(defaultmimzA - defaultmimzB) <= 0.1;
        agreeAC = abs(defaultmimzA - defaultmimzC) <= 0.1;
        agreeBC = abs(defaultmimzB - defaultmimzC) <= 0.1;
        defaultmimz = mz;
        if agreeAB && agreeAC && agreeBC
            defaultmimz = defaultmimzA;
            decision = 'ABC';
        elseif agreeAB
            defaultmimz = defaultmimzA;
            decision = 'AB';
        elseif agreeBC
            defaultmimz = defaultmimzB;
            decision = 'BC';
        elseif agreeAC
            defaultmimz = defaultmimzA;
            decision = 'AC';
        end

        if abs(defaultmimz - mz) < 0.1
            defaultmimz = mz;
        else
            spectrumsectionind = abs(spectrumout(:,1) - defaultmimz) <= 0.05;
            finalpeak = spectrumout(spectrumsectionind,:);
            defaultmimz = sum(finalpeak(:,1).*finalpeak(:,2))/sum(finalpeak(:,2));
        end
    end
end
end