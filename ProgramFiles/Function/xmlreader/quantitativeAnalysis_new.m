function msdata_QA = quantitativeAnalysis_new(scannums,msdata,MS1tol,MS1tolUnit,MS2tol,MS2tolUnit,simscoretol,chnos,areayn,mergeyn,plotyn,distyn)
%QUANTITATIVEANALYSIS: This function performs quantitative analysis on a
% given msdata structure. 
%
% Syntax: 
%   msdata_QA = quantitativeAnalysis(scannums,msdata,MS1tol,MS1tolUnit,MS2tol,MS2tolUnit,simscoretol,chnos,areayn,mergeyn,plotyn,distyn)
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
% SCANCLUSTERS, SCOREAUC, MERGEMSN, THINMSDATA
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
        end
        thisretime = allretime(allscannum == allprescan{scannumind}); % retention time of corresponding MS1 scan
        [~,thisretimems1ind] = min(abs(ms1retime - thisretime)); % find retention time index
        leftrtind = thisretimems1ind; % initial left index for retention range
        rightrtind = thisretimems1ind + 1; % initial right index for retention range
        TICint = [];
        go_on = true;
        while go_on % determining final left index
            leftspec = ms1spectra{leftrtind};
            if distyn == 1 % if using isotopic mass distribution to calculate TIC
                precpkind_dist = (abs(leftspec(:,1) - precursorMz) <= MS1tolC); % ******************************************
                for oo = 1:length(massdist)
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
            else
                go_on = false;
            end
            TICint = [leftint;TICint];
        end
        if rightrtind <= length(ms1spectra)
            go_on = true;
            while go_on % determining final right index
                rightspec = ms1spectra{rightrtind};
                if distyn == 1 % if using isotopic mass distribution to calculate TIC
                    precpkind_dist = (abs(rightspec(:,1) - precursorMz) <= MS1tolC); % ******************************************
                    for pp = 1:length(massdist)
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
                    go_on = false;
                end
                TICint = [TICint;rightint];
            end
        else
            rightrtind = thisretimems1ind;
        end
        select_TIC = [reshape(ms1retime(leftrtind:rightrtind),[],1),TICint]; % selecting area of TIC between left and right index
        retention_range{scannumind} = zeros(1,2); % specifying size of retention range
        % Finding local minima to determine correct isomeric section of TIC
        TF = islocalmin(select_TIC(:,2));               % 2018b feature!!
        splits = nnz(TF);                               % number of non-zero entries
        if splits == 0
            retention_range{scannumind}(1,1) = select_TIC(1,1);
            retention_range{scannumind}(1,2) = select_TIC(end,1);
        elseif splits == 1
            retmin = select_TIC(TF);
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
            scan_clusters(count,2) = retention_range{scannumind}(1,1);
            scan_clusters(count,3) = retention_range{scannumind}(1,2);
            scan_clusters(count,4) = allscannum(scannumind);
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
    scan_clusters = scanClusters(scan_clusters,MS1tolC,allscannum,allmslvl,allretime,allprecursormz); % determine scan clusters
    msdata_QA.scanClusters = scan_clusters;
    sim_score = scoreAUC(scan_clusters,allspectra,MS2tol,MS2tolUnit); % calculate similarity scores for scans in each cluster
    msdata_QA.simScore = sim_score;
    msdata_QA = mergeMSn(msdata_QA,simscoretol); % merge scans in clusters that mass similarity score tolerance
    msdata_QA = thinmsdata(msdata_QA);
end
if plotyn == 1 % if desired, plot similarity scores
    simPlot(scannums,msdata);
end
end
