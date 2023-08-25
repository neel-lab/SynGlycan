function msdata = defaultmonoiso(scannums,msdata,defaultisodist,chnos,distyn,varargin)
%DEFAULTMONOISO: This function facilitates Averagine model calculations to
%     correct the precursor m/z values for all MSn scans based on a supplied
%     isotopic distribution.
%
% Syntax: 
% [msdata] = defaultmonoiso(scannums,msdata,defaultisodist,chnos,distyn)
% [msdata] = defaultmonoiso(scannums,msdata,defaultisodist,chnos,distyn,SASSO)
%
% Input: 
% scannums: n x 1 numerical array. Scan numbers to be analyzed.
% msdata: structure. MS experiment data.
% chnos: 1 x 5 numerical array. User custom averagine model. 
% distyn: logical. If true, in theoretical isotope distribution, only the
%     isotopic species with abundance >= 10% will be used.
% SASSO: m x n numerical array. Scan association in triggered experiment.
%     See BUILDM2SASSO for details.
% 
% Output: 
% msdata: structure. Experiment data with precursor monoisotopic mass
%     adjusted.
%
% Note:
% Input "chnos" must be normzlied to 1 Da: Suppose chnos = [a,b,c,d,e],
%     then [AM_C,AM_H,AM_N,AM_O,AM_S] .* chnos must equals to 1.
%     Here "AM_C" means atomic mass of C.
% If input "SASSO" is provided, to save time the calculation will be carried
%     out on triggering scans only. Because triggered scans share the
%     precursor ion with triggering scans, there is no need to repeat
%     calculation.
% 
% Examples:
% [msdata] = defaultmonoiso(msdata.scannum,msdata,default_iso_dist,[],1);
% 
% Children function: 
% GETDEFAULTMONOISOWDIST
% 
% See also: 
% PREPROCESSGUI  DEFAULTMONOISO  QUANTITATIVEANALYSIS
%

h = waitbar(0,'','Name','Performing Averagine Calculation...');
allscannum = msdata.scannum;
allcharge = msdata.charge;
allretime = msdata.retime;
allmslevel = msdata.mslvl;
allspectra = msdata.spectra;
allmz = msdata.precursormz;
allprescan = msdata.precursorScanNum;
charge = msdata.charge;
defaultmimz = zeros(length(allscannum),1);
mzmassdist = cell(length(allscannum),1);
decision = cell(length(allscannum),1);
dist_range = defaultisodist{2,1}(1,1)-defaultisodist{1,1}(1,1);
if ~isempty(varargin)  % Triggered
    SASSO = varargin{1};
    for ii = 1:length(scannums)  % MS2 scan number
        scannum_idx = find(allscannum == scannums(ii));
        if ismember(scannums(ii),SASSO(:,1))
            [defaultmimz(scannum_idx),mzmassdist{scannum_idx},decision{scannum_idx}] = getdefaultmonoisowdist(...
                allprescan{scannum_idx},allmz(scannum_idx),charge(scannum_idx),...
                allscannum,allcharge,allretime,allmslevel,allspectra,defaultisodist,dist_range,...
                distyn,chnos);
        elseif ismember(scannums(ii),SASSO(:,2:end-2))
            [HCD,~] = find(SASSO(:,2:end-2)==scannums(ii));
            HCDscan_idx = allscannum == SASSO(HCD,1);
            defaultmimz(scannum_idx) = defaultmimz(HCDscan_idx);
            mzmassdist{scannum_idx} = mzmassdist{HCDscan_idx};
        else
            defaultmimz(scannum_idx) = -1;
            mzmassdist{scannum_idx} = -1;
        end
        if mod(ii,200) == 0
            waitbar(ii/length(scannums),h,sprintf('%d/%d',ii,length(scannums)));
        end
    end  
else  % Non-triggered
    for ii = 1:length(scannums)
        scannum_idx = find(allscannum == scannums(ii));
        if allmslevel(scannum_idx) > 1
            [defaultmimz(scannum_idx),mzmassdist{scannum_idx},decision{scannum_idx}] = getdefaultmonoisowdist(...
                allprescan{scannum_idx},allmz(scannum_idx),charge(scannum_idx),...
                allscannum,allcharge,allretime,allmslevel,allspectra,defaultisodist,dist_range,...
                distyn,chnos);
        else
            defaultmimz(scannum_idx) = -1;
            mzmassdist{scannum_idx} = -1;
        end
        if mod(ii,200) == 0
            waitbar(ii/length(scannums),h,sprintf('%d/%d',ii,length(scannums)));
        end
    end    
end 
msdata.originalprecursormz = allmz;
msdata.precursormz = defaultmimz;
msdata.precursormzmethod = decision;
msdata.defaultisodist = defaultisodist;
if distyn == 1
    msdata.theomassdist = mzmassdist;
end
close(h)
end