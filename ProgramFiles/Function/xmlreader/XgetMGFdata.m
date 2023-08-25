function msdata = XgetMGFdata(xmlobj,infos,h)
%XGETMGFDATA: This function reads selected fields available in MGF
% files
%
% Syntax:
%   msdata = XgetMGFdata(xmlobj,infos,h)
%
% Input: 
%   msfilename and location, information to extract, waitbar from
%   getmsdata
% 
% Output: 
%   msdata structure containing info on: i. scannumber; ii. mslvl;
%   iii. retention time; iv. fragmode and charge; precursor m/s and scannum;
%   v. full spectrum, vi. basepeak, vii. number of data points and viii. TIC.
%
% Examples:
%   msdata = getmzxmldata('a filename',{'std','spectra','totIonCurrent'},h = waitbar(0,'','Name','Writing msdata file...'));
%
%See also: GETMSDATA, GETINDEXEDMZMLDATA, GETMZMLDATA, GETMZXMLDATA, PROTEOWIZARD
%

standardinfo = {'num','msLevel','retentionTime',...
    'precursorMz','precursorCharge','fragmode'};
additionalinfo = {'peaksCount','polarity','scanType',...
    'centroided','deisotoped','chargeDeconvoluted',...
    'lowMz','highMz','basePeakMz','basePeakIntensity',...
    'totIonCurrent','precursorScanNum','ionisationEnergy',...
    'precision'}; %'startMz','endMz','spectra'
manufacturer = {'Thermo'};
if nargin == 0
    disp(strjoin(['Standard: ',standardinfo],', '))
    disp(strjoin(['Custom: ',additionalinfo],', '))
    disp(['Default MS inst. manufacturer: ',manufacturer{1}])
    msdata = [];
else
    
    if isempty(infos)
        infos = {'std'};
    elseif ischar(infos)
        infos = {infos};
    end
    stdinfo = standardinfo(ismember(lower(standardinfo),lower(infos)));
    addinfo = additionalinfo(ismember(lower(additionalinfo),lower(infos)));
    
    ion_tag = regexp(xmlobj,'BEGIN IONS');
    tagClose = regexp(xmlobj,'END IONS');
    ms2_count = length(ion_tag);
    last = xmlobj(ion_tag(end):(tagClose(end)+7));
    last_text = extractBetween(last,'.','.');
    scancount = str2double(last_text{1});
    
    if ismember('std',infos)
        scanum = zeros(scancount,1);
        retime = cell(size(scanum));
        precursormz = zeros(size(scanum));
        charge = zeros(scancount,1);
    end
    if ~isempty(stdinfo)
        stdinfosto = cell(scancount,length(stdinfo));
    end
    if ~isempty(addinfo)
        addinfosto = cell(scancount,length(addinfo));
    end

    for ii = 1:length(ion_tag)
        scan_text = xmlobj((ion_tag(ii)):(tagClose(ii)+7));
        idx = str2double(extractBetween(scan_text,'.','.'));
        idx = idx(1);
        if ismember('std',infos)
            try
                scanum(idx) = idx;
                retime(idx) = str2double(extractBetween(scan_text,'RTINSECONDS=',' '));
                precursormz(idx) = str2double(extractAfter(scan_text,'PEPMASS=',' '));
                charge(idx) = str2double(extractBetween(scan_text,'CHARGE=',' '));
            catch
                scanum(idx) = 0;
                retime(idx) = [];
                precursormz(idx) = [];
                charge(idx) = -1;
            end
        end
        if ~isempty(stdinfo)
            for jj = 1:length(stdinfo)
                if strcmpi(stdinfo{jj},'retentionTime')
                    stdinfosto{idx,jj} = str2double(extractBetween(scan_text,'RTINSECONDS=',' '));
                elseif strcmpi(stdinfo{jj},'num')
                    stdinfosto{idx,jj} = idx;
                elseif strcmpi(stdinfo{jj},'precursorMz')
                    stdinfosto{idx,jj} = str2double(extractAfter(scan_text,'PEPMASS=',' '));
                elseif strcmpi(stdinfo{jj},'precursorCharge')
                    stdinfosto{idx,jj} = str2double(extractBetween(scan_text,'CHARGE=',' '));
                end
            end
        end        
        if ~isempty(addinfo)
            for jj = 1:length(addinfo)
%                 elseif strcmpi(addinfo{j},'peaksCount')
%                     addinfosto{idx,j} = ;
%                 elseif strcmpi(addinfo{j},'startMz')
%                     addinfosto{i,j} = ;
%                 elseif strcmpi(addinfo{j},'endMz')
%                     addinfosto{i,j} = ;
%                 elseif strcmpi(addinfo{j},'lowMz')
%                     addinfosto{idx,j} = ;
%                 elseif strcmpi(addinfo{j},'highMz')
%                     addinfosto{idx,j} = ;
%                 elseif strcmpi(addinfo{j},'basePeakMz')
%                     addinfosto{idx,j} = ;
%                 elseif strcmpi(addinfo{j},'basePeakIntensity')
%                     addinfosto{idx,j} = ;
                if strcmpi(addinfo{jj},'spectra')
                    breaks = regexp(xmlobj, '[\n]');
                    eq = regexp(scan_text,'=');
                    leq = eq(end);
                    b = find(breaks>leq,1);
                    spectra = scan_text(breaks(b)+1:tagClose(ii)-1);
                    columns = textscan(spectra,'%.f %f');
                    mz = columns{1,1};
                    int = columns{1,2};
                    addinfosto{idx,jj} = [mz,int];
                end
            end
        end
        if mod(idx,200) == 0
            waitbar(idx/scancount,h,sprintf('%d/%d',idx,scancount));
        end
    end
    if ismember('std',infos)
        msdata.scannum = scanum;
        msdata.retime = retime;
        msdata.precursormz = precursormz;
        msdata.charge = charge;
    end
    if ~isempty(stdinfo)
        for ii = 1:length(addinfo)
            eval(['msdata.',stdinfo{ii},'=stdinfosto(:,',num2str(ii),');']);
        end
    end
    if ~isempty(addinfo)
        for ii = 1:length(addinfo)
            eval(['msdata.',addinfo{ii},'=addinfosto(:,',num2str(ii),');']);
        end
    end
    close(h);
end