function msdata_QA = quantitativeAnalysis(scannums,msdata,MS1tol,MS1tolUnit,MS2tol,simscoretol,px,retlim,chnos,areayn,mergeyn,plotyn,distyn)
%QUANTITATIVEANALYSIS: This function performs quantitative analysis on a
% given msdata structure. 
%
% Syntax: 
%   msdata_QA = quantitativeAnalysis(scannums,msdata,MS1tol,MS1tolUnit,MS2tol,simscoretol,chnos,areayn,mergeyn,plotyn,distyn)
%
% Input: 
%   scan numbers, msdata structure, MS1 mass tolerance and units, MS2
%   mass tolerance and units, similarity score cut-off, CHNOS ratio for
%   isotopic distribution, optional parameters
%
% Output: 
%   new msdata structure containing info on: i. selectTIC, ii.
%   parsedAUC, iii. retentionRanges, iv. simScore, v. scanClusters
%
% Examples:
%   msdata = quantitativeAnalysis(msdata.scannum,msdata,10,'ppm',1,'Da',0.6,[],1,1,0,1);
%
%See also: PREPROCESSGUI, GETMSDATA, MASSDISTRIBUTION, CALCULATEAUC,
% SCANCLUSTERS, SCOREAUC, SPLITSCOMPARE, MERGEMSN, THINMSDATA
%

f = waitbar(0,'','Name','Performing Quantitative Analysis...');
% Variables for all scans from msdata
msdata_QA = msdata;
allmslvl = msdata.mslvl;
allcharge = msdata.charge;
allretime = msdata.retime;
allscannum = msdata.scannum;
allspectra = msdata.spectra;
allprecursormz = msdata.precursormz;
allprescan = msdata.precursorScanNum;
allfragmode = msdata.fragmode;
% Setting mass tolerance if units = Da
if (strcmpi(MS1tolUnit,'Da')) % for Da
    MS1tolC = MS1tol;  % mass tolerance in Da
elseif ~(strcmpi(MS1tolUnit,{'ppm','Da'}))
    error('MATLAB:GLYCOPAT:ERRORUNIT','INCORRECT UNIT FOR MS1 MASS TOL');
end
% Specifying empty variable sizes
retention_range = cell(length(allmslvl),1);
TICdata = cell(length(allmslvl),1);
if areayn == 1 % if area calculations are to be done
    selectAUCdata = zeros(length(allmslvl),1);
end
if mergeyn == 1 % if merging is to be done
    scan_clusters = [0,0,0,0,0];
    count = 1;
end
% Beginning of Calculations
ismslvl1 = allmslvl == 1; % find all ms1 scans
ms1retime = allretime(ismslvl1); % retention times for all ms1 scans
ms1spectra = allspectra(ismslvl1); % spectra for all ms1 scans
for ii = 1:length(scannums)%scannumind = 1:length(allmslvl)
    scannumind = find(allscannum == scannums(ii)); % find index for a specified scan number
    if allmslvl(scannumind) > 1 % for all MSn scans
        precursorMz = allprecursormz(scannumind);
        charge = allcharge(scannumind);
        frag = allfragmode{scannumind};
        % Setting mass tolerance if units = ppm
        if strcmpi(MS1tolUnit,'ppm')  % for ppm
            MS1tolC = MS1tol/1e6*precursorMz;  % mass tolerance in ppm
        end  
        if distyn == 1 % if want to use isotopic mass distribution, obtain m/z values
            if isfield(msdata,'theomassdist')
                massdist = msdata.theomassdist{scannumind}; % values saved from Averagine model
            else
                defaultmonoisodist_mat = load('default_iso_dist.mat'); % load default isotopic distribution
                defaultmonoisodist = defaultmonoisodist_mat.defaultisodist;
                massdist = massDistribution(scannumind,precursorMz,allcharge(scannumind),allretime,allmslvl,allspectra,defaultmonoisodist,chnos); % calculate m/z values
%                 msdata.theomassdist{scannumind} = allmassdist;
            end
        else
            massdist = [];
        end
        thisretime = allretime(allscannum == allprescan{scannumind}); % retention time of corresponding MS1 scan
        [~,thisretimems1ind] = min(abs(ms1retime - thisretime)); % find retention time index
        leftrtind = thisretimems1ind; % initial left index for retention range
        rightrtind = thisretimems1ind + 1; % initial right index for retention range
        TICint = [];
        [o,~] = size(massdist);
        go_on = true;
        addTIC = true;
        bump = 0;
        while go_on % determining final left index
            leftspec = ms1spectra{leftrtind};
            if distyn == 1 % if using isotopic mass distribution to calculate TIC
                precpkind_dist = (abs(leftspec(:,1) - precursorMz) <= MS1tolC); % ******************************************
                for oo = 1:o
                    aa = abs(leftspec(:,1) - massdist(oo,1)) <= MS1tolC;
                    precpkind_dist = [precpkind_dist,aa];
                end
                precpkind = any(precpkind_dist == 1,2);
            else % if using only precursor m/z to calculate TIC
                precpkind = (abs(leftspec(:,1) - precursorMz) <= MS1tolC);
            end
            leftint = sum(leftspec(precpkind,2));
            if leftint > 0 && leftrtind > 1
                leftrtind = leftrtind - 1;
            elseif leftrtind == thisretimems1ind
                TICint = [leftint;TICint];
                if thisretimems1ind < 150
                    lepm = thisretimems1ind-1;
                else
                    lepm = 150;
                end
                for ff = 1:lepm                  
                    leftrtind = thisretimems1ind - ff;
                    leftspec = ms1spectra{leftrtind};
                    if distyn == 1 % if using isotopic mass distribution to calculate TIC
                        precpkind_dist = (abs(leftspec(:,1) - precursorMz) <= MS1tolC); % ******************************************
                        for oo = 1:o
                            aa = abs(leftspec(:,1) - massdist(oo,1)) <= MS1tolC;
                            precpkind_dist = [precpkind_dist,aa];
                        end
                        precpkind = any(precpkind_dist == 1,2);
                    else % if using only precursor m/z to calculate TIC
                        precpkind = (abs(leftspec(:,1) - precursorMz) <= MS1tolC);
                    end
                    leftint = sum(leftspec(precpkind,2));
                    if leftint > 0 && leftrtind > 1
                        leftrtind = leftrtind - 1;
                        break
                    end
                    TICint = [leftint;TICint];
                    if ff == lepm
                        addTIC = false;
                    end
                end
            else
                if bump == 0
                    if leftrtind ~= 1
                        bump = 1;
                        leftrtind = leftrtind - 1;
                    else
                        go_on = false;
                    end
                else
                    go_on = false;
                end
            end
            if addTIC
                TICint = [leftint;TICint];
            end
        end
        if rightrtind <= length(ms1spectra)
            go_on = true;
            addTIC = true;
            bump = 0;
            while go_on % determining final right index
                rightspec = ms1spectra{rightrtind};
                if distyn == 1 % if using isotopic mass distribution to calculate TIC
                    precpkind_dist = (abs(rightspec(:,1) - precursorMz) <= MS1tolC); % ******************************************
                    for pp = 1:o
                        bb = abs(rightspec(:,1) - massdist(pp,1)) <= MS1tolC;
                        precpkind_dist = [precpkind_dist,bb];
                    end
                    precpkind = any(precpkind_dist == 1,2);
                else % if using only precursor m/z to calculate TIC
                    precpkind = (abs(rightspec(:,1) - precursorMz) <= MS1tolC);
                end
                rightint = sum(rightspec(precpkind,2));
                if rightint > 0 && rightrtind < length(ms1retime)
                    rightrtind = rightrtind + 1;
                else
                    if bump == 0
                        if rightrtind ~= length(ms1spectra)
                            rightrtind = rightrtind + 1;
                            bump = 1;
                        else
                            go_on = false;
                        end
                    else
                        go_on = false;
                    end
                end
                if addTIC
                    TICint = [TICint;rightint];
                end
            end
        else
            rightrtind = thisretimems1ind;
        end
        try
        select_TIC = [reshape(ms1retime(leftrtind:rightrtind),[],1),TICint]; % selecting area of TIC between left and right index
        catch
            a=1;
        end
        retention_range{scannumind} = zeros(1,2); % specifying size of retention range
        % Finding local minima to determine correct isomeric section of TIC
        [splits,cuts] = splitsCompare(select_TIC,px);
%         TF = islocalmin(select_TIC(:,2));               % 2018b feature!!
%         splits = nnz(TF);                               % number of non-zero entries
        if size(splits,1) == 0
            retention_range{scannumind}(1,1) = select_TIC(1,1);
            retention_range{scannumind}(1,2) = select_TIC(end,1);
        elseif size(splits,1) == 1
            retmin = cuts;
            if (thisretime >= select_TIC(1,1)) && (thisretime <= retmin)
                retention_range{scannumind}(1,1) = select_TIC(1,1);
                retention_range{scannumind}(1,2) = retmin;
            else
                retention_range{scannumind}(1,1) = retmin;
                retention_range{scannumind}(1,2) = select_TIC(end,1);
            end
        else
            retmin = cuts;
            if (thisretime >= select_TIC(1,1)) && (thisretime <= retmin(1))
                retention_range{scannumind}(1,1) = select_TIC(1,1);
                retention_range{scannumind}(1,2) = retmin(1);
            elseif (thisretime == retmin(end))
                retention_range{scannumind}(1,1) = retmin(end-1);
                retention_range{scannumind}(1,2) = retmin(end);
            elseif (thisretime > retmin(end)) && (thisretime <= select_TIC(end,1))
                retention_range{scannumind}(1,1) = retmin(end);
                retention_range{scannumind}(1,2) = select_TIC(end,1);
            else
                for rr = 2:size(splits,1)
                    if (thisretime >= retmin(rr-1)) && (thisretime < retmin(rr))
                        retention_range{scannumind}(1,1) = retmin(rr-1);
                        retention_range{scannumind}(1,2) = retmin(rr);
                        break
                    end
                end
            end
        end
        if retention_range{scannumind}(1,1) == 0
            error('Retention Range = 0')
        end
        if areayn == 1 % if desired, calculate area under TIC curve
            selectAUCdata(scannumind) = calculateAUC(select_TIC,retention_range{scannumind});
        end
        TICdata{scannumind} = select_TIC;
        if mergeyn == 1 % if merging is desired, save information for further calculations
            scan_clusters(count,1) = precursorMz;
            scan_clusters(count,2) = charge;
            if strcmpi(frag,'CID')
                scan_clusters(count,3) = 1;
            elseif strcmpi(frag,'HCD')
                scan_clusters(count,3) = 2;
            elseif strcmpi(frag,'ETD')
                scan_clusters(count,3) = 3;
            elseif strcmpi(frag,'ETciD')
                scan_clusters(count,3) = 4;
            elseif strcmpi(frag,'EThcD')
                scan_clusters(count,3) = 5;
            else
                scan_clusters(count,3) = 0;
            end
            scan_clusters(count,4) = retention_range{scannumind}(1,1);
            scan_clusters(count,5) = retention_range{scannumind}(1,2);
            scan_clusters(count,6) = allscannum(scannumind);
            count = count + 1;
        end
    end
    waitbar(ii/length(scannums),f,sprintf('%d/%d',ii,length(scannums)));
end
msdata_QA.retentionRange = retention_range;
msdata_QA.selectTIC = TICdata;
if areayn == 1
    msdata_QA.parsedAUC = selectAUCdata;
end
close(f)
if mergeyn == 1 % if desired, perform merging
    % >>> Cluster scans based on precursor m/z and retention time range
    scan_clusters = scanClusters(scan_clusters,MS1tolC,allscannum,allmslvl,allretime,allprecursormz,allcharge,allfragmode); % determine scan clusters
    %msdata_QA.scanClusters = scan_clusters;
    sim_score = scoreAUC(scan_clusters,allscannum,allspectra,MS2tol); % calculate similarity scores for scans in each cluster
    %msdata_QA.simScore = sim_score;
    msdata_QA = mergeMSn1(msdata_QA,sim_score,simscoretol,MS2tol); % merge scans in clusters that mass similarity score tolerance
    msdata_QA = thinmsdata(msdata_QA);
    allmslvl = msdata_QA.mslvl;
    allcharge = msdata_QA.charge;
    allretime = msdata_QA.retime;
    allscannum = msdata_QA.scannum;
    allspectra = msdata_QA.spectra;
    allprecursormz = msdata_QA.precursormz;
    allfragmode = msdata_QA.fragmode;
    allTIC = msdata_QA.selectTIC;
    scan_clusters = [0,0,0,0,0];
    count = 1;
    for xx = 1:length(allmslvl)
        if allmslvl(xx) > 1
            if ~isempty(allTIC{xx})
                scan_clusters(count,1) = allprecursormz(xx);
                scan_clusters(count,2) = allcharge(xx);
                if strcmpi(allfragmode{xx},'CID')
                    scan_clusters(count,3) = 1;
                elseif strcmpi(allfragmode{xx},'HCD')
                    scan_clusters(count,3) = 2;
                elseif strcmpi(allfragmode{xx},'ETD')
                    scan_clusters(count,3) = 3;
                elseif strcmpi(allfragmode{xx},'ETciD')
                    scan_clusters(count,3) = 4;
                elseif strcmpi(allfragmode{xx},'EThcD')
                    scan_clusters(count,3) = 5;
                else
                    scan_clusters(count,3) = 0;
                end
                scan_clusters(count,4) = allTIC{xx}(1,1);
                scan_clusters(count,5) = allTIC{xx}(end,1);
                scan_clusters(count,6) = allscannum(xx);
                count = count + 1;
            end
        end
    end
    scan_clusters = scanClusters(scan_clusters,MS1tolC,allscannum,allmslvl,allretime,allprecursormz,allcharge,allfragmode);
    sim_score = scoreAUC(scan_clusters,allscannum,allspectra,MS2tol);
    msdata_QA = mergeMSn2(msdata_QA,sim_score,simscoretol,MS2tol); % merge scans in clusters that mass similarity score tolerance
    msdata_QA = thinmsdata(msdata_QA);
end
if plotyn == 1 % if desired, plot similarity scores
    simPlot(scannums,msdata);
end
end

% a = [scan_clusters(2,:);scan_clusters(find(abs(scan_clusters(:,1) - scan_clusters(2,1)) <= 1),:)];
% b = [scan_clusters(2,:);scan_clusters(find(a(:,2) == scan_clusters(2,2)),:)];
% c = [scan_clusters(2,:);scan_clusters(find(b(:,3) == scan_clusters(2,3)),:)];

function theomassdist = massDistribution(scannumind,mz,chg,allretime,allmslevel,allspectra,defaultmonoisodist,chnos)
%MASSDISTRIBUTION: This function calculates the theoretical m/z mass 
% distribution for a select scan number using the Averagine model
%
% Syntax:
%   theomassdist = massDistribution(scannumind,mz,chg,allretime,allmslevel,allspectra,defaultmonoisodist,chnos)
%
% Input: 
%   Scan number index of the desired scan, select m/z, select charge,
%   retention time for all scans, MS level for all scans, spectra for all
%   scans, default monoisotopic distribution, CHNOS ratio values
% 
% Output: 
%   n x 2 Double containing the m/z values and the relative abundance of
%   these values
%
% Examples:
%   theomassdist = massDistribution(651,714.2861,msdata.charge(651),
%   msdata.retime,msdata.mslvl,msdata.spectra,defaultmonoisodist,[])
%
%See also: QUANTITATIVEANALYSIS_SPLIT
%

dist_range = defaultmonoisodist{2,1}(1,1)-defaultmonoisodist{1,1}(1,1);
mass = mz * chg - 1.0078246*chg;
tempretime = allretime(scannumind);
temprange = (abs(allretime-tempretime) <= .25) & allmslevel == 1;
spectrain=allspectra(temprange);
spectrasize = cellfun(@size,spectrain,num2cell(ones(size(spectrain)))); % calculate size of each of these spectra
spectra = zeros(sum(spectrasize),2); % create empty element with size equal to the sum of the elements in all spectra
ind = 1;
for i = 1:length(spectrasize)        % merge all input spectra in a single variable.
    spectra(ind:ind + spectrasize(i) - 1,:) = double(spectrain{i});
    ind = ind + spectrasize(i);
end
spectrumout = spectra(spectra(:,2) > 0,:);
spectrumsectionind = ((spectrumout(:,1) - mz) <= 5/chg) & ((spectrumout(:,1) - mz) >= -3/chg);
spectrumsection = spectrumout(spectrumsectionind,:);  % consider chg >= 2, 3 Da off monoiso max
idx = floor(mass);
if idx == 0 % if the precursor m/z is 0
    error('No precursor m/z and/or charge value present in msdata.')
elseif idx <= defaultmonoisodist{end,1}(1,1) % if the precursor m/z is within the theoretical distribution
    idx_alt = round(idx/dist_range); % find index of mass based on distribution range
    if idx_alt == 0 % if rounding results in 0
        idx_alt = 1; % set value to 1
    elseif idx_alt > length(defaultmonoisodist) % if rounding results in value greater than length of theoretical distribution
        idx_alt = length(defaultmonoisodist); % set to max
    end
    theoisodist = defaultmonoisodist{idx_alt,1}; % select theoretical distribution
else
    % [0.024321959339170,0.476559048284606,0.001461698046926,0.012941719677186,0] GLYCAN CHNOS VALUES
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
expandtheoiso = [expandtheoiso(:),zeros(length(expandtheoiso),1)];
expandtheoiso = [expandtheoiso;theoisodist];
etisotail = theoisodist(end,1)+theoisostep:theoisostep:max(spectrumsection(:,1))+1;
expandtheoiso = [expandtheoiso;[etisotail(:),zeros(length(etisotail),1)]];
expandtheoiso = unique(expandtheoiso,'rows');
% Collecting all peaks in expandtheoiso about 10%
totIC = sum(expandtheoiso(:,2));
theomassdist = [expandtheoiso(:,1),(expandtheoiso(:,2)/totIC)];
intidx = any(theomassdist(:,2) < 0.1,2);
theomassdist(intidx,:) = [];
end

function [splits,cuts] = splitsCompare(select_TIC,px)
%SPLITSCOMPARE: This function performs quantitative analysis on a
% given msdata structure. 
%
% Syntax: 
%   [splits,cuts] = splitsCompare(select_TIC,px)
%
% Input: 
%   Select TIC data between a given retention time range for a desired
%   precursor m/z and a percent to use as a cut-off.
%
% Output: 
%   splits: Intensity values for the local minima that comply with the
%   user-inputted cut-off
%   cuts: Retention times for the local minima that comply with the
%   user-inputted cut-ff.
%
% Examples:
%   [splits,cuts] = splitsCompare(select_TIC,0.2)
%
%See also: PREPROCESSGUI, GETMSDATA, MASSDISTRIBUTION, CALCULATEAUC,
% QUANTITATIVEANALYSIS, SCANCLUSTERS, SCOREAUC, MERGEMSN, THINMSDATA
%

ret = select_TIC(:,1);
int = select_TIC(:,2);
TG = islocalmax(int);
TF = islocalmin(int);
dips = nnz(TF);
a = int(TG);
b = int(TF);
c = ret(TF);
splits = [];
cuts = [];
hght = zeros(dips,1);
ab = zeros(dips,1);
for ii = 1:length(b)
    if ii < length(a)
        hght(ii) = min(a(ii),a(ii+1));
    else
        hght(ii) = a(ii);
    end
    ab(ii) = b(ii)/hght(ii);
    if ab(ii) < (1-px)
        splits = [splits;b(ii)];
        cuts = [cuts;c(ii)];
    end
end
end

function selectAUCdata = calculateAUC(select_TIC,retention_range) % calculate area between left and right index under TIC curve
%SELECTAUCDATA: This function calculates the area between the left and
% right index under the TIC curve for a specified m/z value.
%
% Syntax:
%   selectAUCdata = calculateAUC(select_TIC,retention_range)
%
% Input: 
%   n x 2 double containing the TIC data for the select m/z, right and left
%   index of the retention time for the desired curve
% 
% Output: 
%   Double containing the area value
%
% Examples:
%   selectAUCdata = calculateAUC([10.7105316666670,0;10.7271816666670,
% 201258.947265625;10.7437983333330,279880.815429688;10.7604533333330,
% 617913.988281250;10.7770750000000,478313.134765625;10.7894400000000,
% 271338.918945313;10.8096450000000,88903.5136718750;10.8262583333330,0],
% {[10.7105316666670,10.8262583333330]})
%
%See also: QUANTITATIVEANALYSIS_SPLIT
%

if length(select_TIC(:,1)) > 1
    leftind = find(select_TIC(:,1) == retention_range(1,1));
    rightind = find(select_TIC(:,1) == retention_range(1,2));
    selectAUCdata = trapz(select_TIC(leftind:rightind,1),select_TIC(leftind:rightind,2)); % calculate area using trapezoidal method
else
    selectAUCdata = -1;
end
end

function scan_clusters = scanClusters(scan_clusters,MS1tolC,allscannum,allmslvl,allretime,allprecursormz,allcharge,allfragmode)
%SCANCLUSTERS: This function generates scan clusters for scans that have
%the same precursor m/z and are in the same retention time range.
%
% Syntax:
%   scan_clusters = scanClusters(scan_clusters,MS1tolC,allscannum,allmslvl,allretime,allprecursormz)
%
% Input: 
%   n x 5 double containing the scan cluster information, MS1 tolerance
%   value, scan numbers for all scans, MS levels for all scans, retention
%   times for all scans, and all precursor m/z values
% 
% Output: 
%   n x m double containing the scan cluster information include the scan
%   numbers for all matching spectra
%
% Examples:
%   scan_clusters = scanClusters([695.341674804688,10.7105316666670,
% 10.8262583333330,649,0;714.560124686258,10.7894400000000,10.809645000000,
% 651,0;695.344970703125,10.7105316666670,10.8262583333330,652,0],0.0191,
% msdata.scannum,msdata.mslvl,msdata.retime,msdata.precursormz)
%
%See also: QUANTITATIVEANALYSIS_SPLIT, MERGEMS, SCOREAUC, THINMSDATA
%
%SCANCLUSTERS: This function generates scan clusters for scans that have
%the same precursor m/z and are in the same retention time range.
%
% Syntax:
%   scan_clusters = scanClusters(scan_clusters,MS1tolC,allscannum,allmslvl,allretime,allprecursormz)
%
% Input: 
%   n x 5 double containing the scan cluster information, MS1 tolerance
%   value, scan numbers for all scans, MS levels for all scans, retention
%   times for all scans, and all precursor m/z values
% 
% Output: 
%   n x m double containing the scan cluster information include the scan
%   numbers for all matching spectra
%
% Examples:
%   scan_clusters = scanClusters([695.341674804688,10.7105316666670,
% 10.8262583333330,649,0;714.560124686258,10.7894400000000,10.809645000000,
% 651,0;695.344970703125,10.7105316666670,10.8262583333330,652,0],0.0191,
% msdata.scannum,msdata.mslvl,msdata.retime,msdata.precursormz)
%
%See also: QUANTITATIVEANALYSIS_SPLIT, MERGEMS, SCOREAUC, THINMSDATA
%

h = waitbar(0,'','Name','Generating Spectra Clusters...');
for ww = 1:size(allscannum)
    scannumind = ww;%find(allscannum == scannums(ww));
    if allmslvl(scannumind)>1 % for all MSn scans
        retime = allretime(scannumind);
        precursorMz = allprecursormz(scannumind);
        charge = allcharge(scannumind);
        frag = allfragmode{scannumind};
        if strcmpi(frag,'CID')
            frag_idx = 1;
        elseif strcmpi(frag,'HCD')
            frag_idx = 2;
        elseif strcmpi(frag,'ETD')
            frag_idx = 3;
        elseif strcmpi(frag,'ETciD')
            frag_idx = 4;
        elseif strcmpi(frag,'EThcD')
            frag_idx = 5;
        end
        origparse = scan_clusters(:,1);
        difference = abs(origparse - precursorMz) <= MS1tolC; % look if precursor m/z is already in matrix within mass tolerance
        if any(difference) % if precursor m/z is already in matrix
            kk = find(difference); % row of precursor m/z if already exists
            for jj = 1:length(kk) % for each row that precursor m/z is in
                if allscannum(scannumind) == scan_clusters(kk(jj),6) % if scan number is already present, skip
                    continue
                elseif (retime >= scan_clusters(kk(jj),4)) && (retime <= scan_clusters(kk(jj),5)) % if precursor m/z is within retention range
                    if charge == scan_clusters(kk(jj),2) % if charge is the same as the index scan
                        if frag_idx == scan_clusters(kk(jj),3) % if fragmode is the same as the index scan     
                            [r,c] = size(scan_clusters);
                            massLocation = zeros(r,1);
                            scan_clusters = [scan_clusters, massLocation]; % add column of zeros to matrix
                            d = c + 1;
                            for ll = 1:d
                                if scan_clusters(kk(jj),ll) == 0
                                    scan_clusters(kk(jj),ll) = allscannum(ww); % replace first zero value in row with scan number
                                break
                                end
                            end
                        end
                    end
                end
            end          
        end
    end
    waitbar(ww/length(allscannum),h,sprintf('%d/%d',ww,length(allscannum)));
end
scan_clusters(~any(scan_clusters,2),:) = []; % remove all rows with only zero values
scan_clusters(:,~any(scan_clusters,1)) = []; % remove all columns with only zero values
close(h)
end

function sim_score = scoreAUC(scan_clusters,allscannum,allspectra,MS2tolmat) % scoring scan clusters
%SCOREAUC: This function scores all of the possible scan clusters generated
%using SCANCLUSTERS.
%
% Syntax:
%   sim_score = scoreAUC(scan_clusters,allspectra,MS2tol,MS2tolUnit)
%
% Input: 
%   n x m double containing the scan cluster information include the scan
%   numbers for all matching spectra, spectra for all scans, MS2 tolerance 
%   value along with the units
% 
% Output: 
%   n x 2 cell containing the original scan and their corresponding scan
%   matches with the score for how well these matches compared to the
%   original scan
%
% Examples:
%   sim_score = scoreAUC([714.560124686258,10.7894400000000,
%   10.8096450000000,651,0,0,0;695.344970703125,10.7105316666670,
%   10.8262583333330,652,649,0,0;674.300903320313,10.8762650000000,
%   11.0129716666670,662,0,0,0],msdata.spectra,1,'Da')
%
%See also: QUANTITATIVEANALYSIS_SPLIT, CALCULTEAUC, MERGEMS,
% SCANCLUSTERS, MASSDISTRIBUTION, THINMSDATA
%

g = waitbar(0,'','Name','Performing Similarity Scoring...');
aa = unique(scan_clusters(:,6));
if aa(1) == 0
    aa = aa(2:end);
end
sim_score = cell(length(aa),4);
for mm = 1:length(aa)
    scan_idx = (scan_clusters(:,6) == aa(mm));
    ri = scan_clusters(scan_idx,:);
    lvl = 1;
    row = nonzeros(ri);
    if length(row) > 6
        scans = row(6:end);
        scannum = row(6);
        MS2tol = MS2tolmat{row(3),1};
        MS2tolUnit = MS2tolmat{row(3),2};
        if (strcmpi(MS2tolUnit,'Da')) % for Da
            MS2tolC = MS2tol;  % mass tolerance in Da
        elseif strcmpi(MS2tolUnit,'ppm')  % for ppm
            MS2tolC = MS2tol/1e6*ri(1);  % mass tolerance in ppm
        elseif ~(strcmpi(MS2tolUnit,{'ppm','Da'}))
            error('MATLAB:GLYCOPAT:ERRORUNIT','INCORRECT UNIT FOR MS2 MASS TOL');
        end
        lscn = length(scans);
        sim_score{mm,1} = scannum;
        sim_score{mm,3} = row(3);
        sim_score{mm,4} = ri(1);
        scannum_idx = allscannum == scannum;
        raw_spec1 = ms2centroid(allspectra{scannum_idx},MS2tol,MS2tolUnit,ri(1)); % centroid spectrum
        if isempty(raw_spec1)
            continue
        end
        [~,idx1] = sort(raw_spec1(:,2),'descend');
        spec1 = raw_spec1(idx1,:);
        m1 = padarray(spec1(:,1),30,'post');
        i1 = padarray(spec1(:,2),30,'post');
        sim_score{mm,2} = zeros(lscn-1,2);
        for nn = 2:lscn
            scans_idx = allscannum == scans(nn);
            raw_spec2 = ms2centroid(allspectra{scans_idx},MS2tol,MS2tolUnit,ri(1));
            if isempty(raw_spec2)
                continue
            end
            % Similarity Score Calculation
            [~,idx2] = sort(raw_spec2(:,2),'descend');
            spec2 = raw_spec2(idx2,:);
            m2 = padarray(spec2(:,1),30,'post');
            i2 = padarray(spec2(:,2),30,'post');
            for k = 1:30
                hp(k) = i1(k); % part numerator
                hp_m = m1(k);
                mat = abs(m2-hp_m) <= MS2tolC;
                if any(mat)
                    loc = find(mat);
                    if length(loc) > 1
                        [~,loc] = min(abs(m2-hp_m));
                    end
                    hpi(k) = i2(loc); % part numerator
                else
                    hpi(k) = 0; % part numerator
                end
                numer(k) = hp(k)*hpi(k); % total numerator
                hp_tot(k) = (hp(k))^2; % part denominator
                hq(k) = (i2(k))^2; % part denominator
            end
            sim_score{mm,2}(lvl,1) = scans(nn);
            sim_score{mm,2}(lvl,2) = sum(numer)/sqrt(sum(hp_tot)*sum(hq)); % similarity score calculation
            lvl = lvl+1;
        end
    end
    waitbar(mm/length(aa),g,sprintf('%d/%d',mm,length(aa)));
end
sim_score(all(cellfun('isempty',sim_score),2),:) = []; % remove all empty rows
close(g)
end

function [msdata] = mergeMSn1(msdata,sim_score,simscoretol,MS2tolmat) % merge all MSn spectra that in each cluster that are similar
%MERGEMSN: This function generates a new msdata structure that
%combines and merges scans that are similar
%
% Syntax:
%   [msdata] = mergeMSn(msdata,simscoretol)
%
% Input: 
%   msdata structure and the tolerance limit for the similarity score
% 
% Output: 
%   A new msdata structure where the matched spectra are present in only
%   one entry
%
% Examples:
%   [msdata] = mergeMSn(msdata,0.6)
%
%See also: QUANTITATIVEANALYSIS_SPLIT, SCANCLUSTERS, SCOREAUC, THINMSDATA
%

j = waitbar(0,'','Name','Performing Primary Merging...');
allscannum = msdata.scannum;
consensus = msdata.spectra;
scan_indexes = cell(length(msdata.spectra),1);
for ii = 1 : length(sim_score)
    scan_idx = sim_score{ii,1};
    if scan_idx == 0
        continue
    end
    scan_indexes{allscannum == scan_idx} = scan_idx;
    sim_score1 = sim_score{ii,2};
    if ~isempty(sim_score1)
%         if size(sim_score1,1) > 1
%             sim_std = std(sim_score1(:,2));
%             if sim_std <= simscoretol*0.1 % Use standard deviation to accept or reject entire groups of scans with close similarity scores
%                 simscoretol_mean = mean(sim_score1(:,2));
%                 sim_score1(:,2) = simscoretol_mean;
%             end
%         end
        for jj = 1:size(sim_score1(:,1),1)
            if sim_score1(jj,2) >= simscoretol
                if isempty(consensus{allscannum == sim_score1(jj,1)})
                    continue
                end
                combine_spec = [consensus{allscannum == scan_idx};msdata.spectra{allscannum == sim_score1(jj,1)}]; % add spectra in a cluster
                [~,idx] = sort(combine_spec(:,1)); % sort combined spectrum based on m/z
                consensus{allscannum == scan_idx} = combine_spec(idx,:);
                scan_indexes{allscannum == scan_idx} = [scan_indexes{allscannum == scan_idx};sim_score1(jj,1)];
                consensus{allscannum == sim_score1(jj,1)} = []; % removing spectra for sim_score1(jj,1)
                score_idx = cell2mat(sim_score(:,1)) == sim_score1(jj,1);
                if sum(score_idx) ~= 0
                    sim_score{score_idx,1} = 0;
                end
            end
        end
        MS2tol = MS2tolmat{sim_score{ii,3},1};
        MS2tolUnit = MS2tolmat{sim_score{ii,3},2};
        consensus{allscannum == scan_idx} = ms2centroid(consensus{allscannum == scan_idx},MS2tol,MS2tolUnit,sim_score{ii,4}); % centroid spectrum
    end
    waitbar(ii/length(sim_score),j,sprintf('%d/%d',ii,length(sim_score)));
end
msdata.indexes = scan_indexes;
msdata.spectra = consensus;
close(j)
end

function [msdata] = mergeMSn2(msdata,sim_score,simscoretol,MS2tolmat) % merge all MSn spectra that in each cluster that are similar
%MERGEMSN: This function generates a new msdata structure that
%combines and merges scans that are similar
%
% Syntax:
%   [msdata] = mergeMSn(msdata,simscoretol)
%
% Input: 
%   msdata structure and the tolerance limit for the similarity score
% 
% Output: 
%   A new msdata structure where the matched spectra are present in only
%   one entry
%
% Examples:
%   [msdata] = mergeMSn(msdata,0.6)
%
%See also: QUANTITATIVEANALYSIS_SPLIT, SCANCLUSTERS, SCOREAUC, THINMSDATA
%

j = waitbar(0,'','Name','Performing Secondary Merging...');
allscannum = msdata.scannum;
consensus = msdata.spectra;
scan_indexes = msdata.indexes;
selectTIC = msdata.selectTIC;
retrange = msdata.retentionRange;
for ii = 1 : length(sim_score)
    scan_idx = sim_score{ii,1};
    if scan_idx == 0
        continue
    end
    scan_indexes{allscannum == scan_idx} = [scan_indexes{allscannum == scan_idx}; scan_idx];
    sim_score1 = sim_score{ii,2};
    if ~isempty(sim_score1)
%         if size(sim_score1,1) > 1
%             sim_std = std(sim_score1(:,2));
%             if sim_std <= simscoretol*0.1 % Use standard deviation to accept or reject entire groups of scans with close similarity scores
%                 simscoretol_mean = mean(sim_score1(:,2));
%                 sim_score1(:,2) = simscoretol_mean;
%             end
%         end
        for jj = 1:size(sim_score1(:,1),1)
            if sim_score1(jj,2) >= simscoretol
                if isempty(consensus{allscannum == sim_score1(jj,1)})
                    continue
                end
                combine_spec = [consensus{allscannum == scan_idx};msdata.spectra{allscannum == sim_score1(jj,1)}]; % add spectra in a cluster
                [~,idx] = sort(combine_spec(:,1)); % sort combined spectrum based on m/z
                consensus{allscannum == scan_idx} = combine_spec(idx,:);
                if size(scan_indexes{allscannum == sim_score1(jj,1)},1) > 1
                    scan_indexes{allscannum == scan_idx} = [scan_indexes{allscannum == scan_idx};scan_indexes{allscannum == sim_score1(jj,1)}];
                else
                    scan_indexes{allscannum == scan_idx} = [scan_indexes{allscannum == scan_idx};sim_score1(jj,1)];
                end
                consensus{allscannum == sim_score1(jj,1)} = []; % removing spectra for sim_score1(jj,1)
                score_idx = cell2mat(sim_score(:,1)) == sim_score1(jj,1);
                if sum(score_idx) ~= 0
                    sim_score{score_idx,1} = 0;
                end
            end
        end
        MS2tol = MS2tolmat{sim_score{ii,3},1};
        MS2tolUnit = MS2tolmat{sim_score{ii,3},2};
        consensus{allscannum == scan_idx} = ms2centroid(consensus{allscannum == scan_idx},MS2tol,MS2tolUnit,sim_score{ii,4}); % centroid spectrum
        retrange{allscannum == scan_idx}(1,1) = selectTIC{allscannum == scan_idx}(1,1);
        retrange{allscannum == scan_idx}(1,2) = selectTIC{allscannum == scan_idx}(end,1);
    end
    scan_indexes{allscannum == scan_idx} = unique(scan_indexes{allscannum == scan_idx});
    waitbar(ii/length(sim_score),j,sprintf('%d/%d',ii,length(sim_score)));
end
msdata.indexes = scan_indexes;
msdata.spectra = consensus;
msdata.retentionRange = retrange;
close(j)
end

function msdata = thinmsdata(msdata)
%THINMSDATA: This function removes redundant information from the msdata
%structure for merged scans.
%
% Syntax:
%   msdata = thinmsdata(msdata)
%
% Input: 
%   msdata structure
% 
% Output: 
%   A new msdata structure that removes empty msdata fields that were
%   blanked during MERGEMS
%
% Examples:
%   msdata = thinmsdata(msdata)
%
%See also: QUANTITATIVEANALYSIS_SPLIT, MERGEMS, SCOREAUC, SCANCLUSTERS
%

u = waitbar(0,'','Name','Thinning msdata...');
fields = fieldnames(msdata);
len = length(msdata.spectra);
for ll = 0:(len-1)
    if isempty(msdata.spectra{len-ll})
        hold_len = length(msdata.scannum);
        for ii = 1:length(fields)
            if length(msdata.(fields{ii})) == hold_len
                msdata.(fields{ii})(len-ll) = [];
            end
        end
    end
    waitbar((ll+1)/len,u,sprintf('%d/%d',(ll+1),len));
end
close(u)
end
