function msdata = getIndexedmzMLdata(xmlobj,infos,h)
%GETINDEXEDMZMLDATA: This function reads selected fields available in
% indexed mzML files
%
% Syntax:
%   msdata = getIndexedmzMLdata(xmlobj,infos,h)
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
%   msdata = getIndexedmzMLdata('a filename',{'std','spectra','totIonCurrent'},h = waitbar(0,'','Name','Writing msdata file...'));
%
% Children function:
% BASE64DECODE
%
% See also:
% GETMSDATA  GETMZMLDATA  GETMZXMLDATA  PROTEOWIZARD
%

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 11/10/2020

standardinfo = {'num','msLevel','retentionTime',...
    'precursorMz','allPrecursorMz','precursorCharge','fragmode'};
additionalinfo = {'peaksCount','polarity',...
    'centroided','deisotoped','chargeDeconvoluted',...
    'startMz','endMz','lowMz',...
    'highMz','basePeakMz','basePeakIntensity',...
    'totIonCurrent','precursorScanNum','ionisationEnergy',...
    'precision','spectra','filterLine'}; % 'scanType'
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
    manuinfoind = ismember(lower(manufacturer),lower(infos));
    if ~any(manuinfoind)
        manuinfo = lower(manufacturer{1});
    else
        manuinfo = manufacturer{manuinfoind};
    end
    scan_tag = regexp(xmlobj,'<spectrum index=');
    tagClose = regexp(xmlobj,'</spectrum');
    scancount = length(scan_tag);
    if ismember('std',infos)
        scanum = zeros(scancount,1);
        mslvl = zeros(size(scanum));
        retime = zeros(size(scanum));
        fragmode = cell(size(scanum));
        precursormz = zeros(size(scanum));
        allprecursormz = cell(size(scanum));
        charge = zeros(size(scanum));
    end
    if ~isempty(stdinfo)
        stdinfosto = cell(scancount,length(stdinfo));
    end
    if ~isempty(addinfo)
        addinfosto = cell(scancount,length(addinfo));
    end
    parent_tag = regexp(xmlobj,'<sourceFile id="');
    parent_close = regexp(xmlobj,'</sourceFile');
    parent_length = length(parent_tag);
    sourceFile = cell(parent_length,1);
    for p = 1:parent_length
        parent_text = xmlobj((parent_tag(p)):(parent_close(p)-1));
        sourceFile{p} = extractBetween(parent_text,'name="','"');
    end
    msdata.sourceFileList = sourceFile;
    for ii = 1:scancount
        spectrum_text = xmlobj((scan_tag(ii)):(tagClose(ii)-1));
        scan_text = extractBetween(spectrum_text,'<spectrum ','</scanList');
        precursor_text = extractBetween(spectrum_text,'<precursorList ','</precursorList');
        peaks_text = extractBetween(spectrum_text,'<binaryDataArrayList ','</binaryDataArrayList');
        if ismember('std',infos)
            scanum(ii) = str2double(extractBetween(scan_text,'scan=','"'));
            mslvl(ii) = str2double(extractBetween(scan_text,'"ms level" value="','"'));
            retime(ii) = str2double(extractBetween(scan_text,'"scan start time" value="','"'));
            if ~isempty(precursor_text)
                allprecursormz{ii}(1,1) = str2double(extractBetween(precursor_text,'"selected ion m/z" value="','"'));
                allprecursormz{ii}(1,2) = str2double(extractBetween(precursor_text,'"isolation window target m/z" value="','"'));
                allprecursormz{ii}(1,3) = str2double(extractBetween(scan_text,'Monoisotopic M/Z:" value="','"'));
                try
                    charge(ii) = str2double(extractBetween(spectrum_text,'"charge state" value="','"'));
                catch
                    charge(ii) = 0;
                end
                if allprecursormz{ii}(1,3) == 0
                    precursormz(ii) = allprecursormz{ii}(1,1);
                else
                    precursormz(ii) = allprecursormz{ii}(1,3);
                end
                thisfragmode = extractBetween(precursor_text,'<activation>','</activation>');
                names = regexp(thisfragmode,' name="');
                dissociations = regexp(thisfragmode,' dissociation');
                frag_count = length(dissociations{1});
                this = thisfragmode{1};
                if frag_count == 1
                    idx1 = names{1}+7;
                    idx2 = dissociations{1};
                    frag_text = this(idx1:idx2);
                    if strcmpi(frag_text,'collision-induced ')
                        fragmode{ii} = 'CID';
                    elseif strcmpi(frag_text,'beam-type collision-induced ')
                        fragmode{ii} = 'HCD';
                    elseif strcmpi(frag_text,'electron transfer ')
                        fragmode{ii} = 'ETD';
                    end
                else
                    for f = 1:frag_count
                        idx1 = (names{1}(f))+7;
                        idx2 = dissociations{1}(f);
                        frag_text = this(idx1:idx2);
                        if contains(frag_text,'supplemental')
                            if strcmpi(frag_text,'supplemental collision-induced ')
                                fragmode{ii} = 'ETciD';
                            elseif strcmpi(frag_text,'supplemental beam-type collision-induced ')
                                fragmode{ii} = 'EThcD';
                            end
                        end
                    end
                end
            else
                precursormz(ii) = -2;
                fragmode{ii} = '';
            end
        end
        if ~isempty(stdinfo)
            for jj = 1:length(stdinfo)
                if strcmpi(stdinfo{jj},'fragmode')
                    thisfragmode = extractBetween(precursor_text,'<activation>','</activation>');
                    names = regexp(thisfragmode,' name="');
                    dissociations = regexp(thisfragmode,' dissociation');
                    frag_count = length(dissociations{1});
                    this = thisfragmode{1};
                    if frag_count == 1
                        idx1 = names{1}+7;
                        idx2 = dissociations{1};
                        frag_text = this(idx1:idx2);
                        if strcmpi(frag_text,'collision-induced ')
                            fragmode{ii} = 'CID';
                        elseif strcmpi(frag_text,'beam-type collision-induced ')
                            fragmode{ii} = 'HCD';
                        elseif strcmpi(frag_text,'electron transfer ')
                            fragmode{ii} = 'ETD';
                        end
                    else
                        for f = 1:frag_count
                            idx1 = (names{1}(f))+7;
                            idx2 = dissociations{1}(f);
                            frag_text = this(idx1:idx2);
                            if contains(frag_text,'supplemental')
                                if strcmpi(frag_text,'supplemental collision-induced ')
                                    fragmode{ii} = 'ETciD';
                                elseif strcmpi(frag_text,'supplemental beam-type collision-induced ')
                                    fragmode{ii} = 'EThcD';
                                end
                            end
                        end
                    end
                    stdinfosto{ii,jj} = fragmode{ii};
                elseif strcmpi(stdinfo{jj},'retentionTime')
                    stdinfosto{ii,jj} = str2double(extractBetween(scan_text,'"scan start time" value="','"'));
                elseif strcmpi(stdinfo{jj},'num')
                    stdinfosto{ii,jj} = str2double(extractBetween(scan_text,'scan=','"'));
                elseif strcmpi(stdinfo{jj},'msLevel')
                    stdinfosto{ii,jj} = str2double(str2double(extractBetween(scan_text,'"ms level" value="','"')));
                elseif strcmpi(stdinfo{jj},'PrecursorMz')
                    try
                        precursormz = str2double(extractBetween(scan_text,'Monoisotopic M/Z:" value="','"'));
                    catch
                        precursormz = str2double(extractBetween(precursor_text,'"selected ion m/z" value="','"'));
                    end
                    if precursormz == 0
                        stdinfosto{ii,jj} = -2;
                    else
                        stdinfosto{ii,jj} = precursormz;
                    end
                elseif strcmpi(stdinfo{jj},'allPrecursorMz')
                    try
                        allprecursormz(1,1) = str2double(extractBetween(precursor_text,'"selected ion m/z" value="','"'));
                        allprecursormz(1,2) = str2double(extractBetween(precursor_text,'"isolation window target m/z" value="','"'));
                        allprecursormz(1,3) = str2double(extractBetween(scan_text,'Monoisotopic M/Z:" value="','"'));
                        stdinfosto{ii,jj} = allprecursormz;
                    catch
                        stdinfosto{ii,jj} = [];
                    end
                elseif strcmpi(stdinfo{jj},'precursorCharge')
                    try
                        stdinfosto{ii,jj} = str2double(extractBetween(precursor_text,'"charge state" value="','"'));
                    catch
                        stdinfosto{ii,jj} = [];
                    end
                end
            end
        end
        if ~isempty(addinfo)
            for jj = 1:length(addinfo)
                addinfosto{ii,jj} = -1;
                if strcmpi(addinfo{jj},'peaksCount')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'defaultArrayLength="','"'));
                elseif strcmpi(addinfo{jj},'polarity')
                    if ismember(scan_text{1},'positive scan')
                        addinfosto{ii,jj} = '+';
                    elseif ismember(scan_text{1},'negative scan')
                        addinfosto{ii,jj} = '-';
                    else
                        addinfosto{ii,jj} = [];
                    end
                elseif strcmpi(addinfo{jj},'centroided')
                    if ismember(scan_text{1},'centroid spectrum')
                        addinfosto{ii,jj} = 'True';
                    else
                        addinfosto{ii,jj} = [];
                    end
                elseif strcmpi(addinfo{jj},'deisotoped')
                    if ismember(scan_text{1},'deisotoped')
                        addinfosto{ii,jj} = 'Yes';
                    else
                        addinfosto{ii,jj} = [];
                    end
                elseif strcmpi(addinfo{jj},'chargeDeconvoluted')
                    if ismember(scan_text{1},'chargeDeconvoluted')
                        addinfosto{ii,jj} = 'Yes';
                    else
                        addinfosto{ii,jj} = [];
                    end
                elseif strcmpi(addinfo{jj},'startMz')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"scan window lower limit" value="','"'));
                elseif strcmpi(addinfo{jj},'endMz')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"scan window upper limit" value="','"'));
                elseif strcmpi(addinfo{jj},'lowMz')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"lowest observed m/z" value="','"'));
                elseif strcmpi(addinfo{jj},'highMz')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"highest observed m/z" value="','"'));
                elseif strcmpi(addinfo{jj},'basePeakMz')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"base peak m/z" value="','"'));
                elseif strcmpi(addinfo{jj},'basePeakIntensity')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"base peak intensity" value="','"'));
                elseif strcmpi(addinfo{jj},'totIonCurrent')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"total ion current" value="','"'));
                elseif strcmpi(addinfo{jj},'precursorScanNum')
                    if any(strfind(spectrum_text,'master scan number'))
                        temp = str2double(extractBetween(spectrum_text,'name="master scan number" value="','"'));
                        if ~isempty(temp)
                            addinfosto{ii,jj} = temp;
                        end
                    else
                        try
                            temp = str2double(extractBetween(precursor_text,'scan=','"'));
                            if ~isempty(temp)
                                addinfosto{ii,jj} = temp;
                            end
                        catch
                            addinfosto{ii,jj} = -1;
                        end
                    end
                elseif strcmpi(addinfo{jj},'ionisationEnergy')
                    try
                        thisfragmode = extractBetween(precursor_text,'<activation>','</activation>');
                        names = regexp(thisfragmode,' name="');
                        dissociations = regexp(thisfragmode,' dissociation');
                        d = length(dissociations{1});
                        values = regexp(thisfragmode,'" value="');
                        refs = regexp(thisfragmode,'" unitCvRef="');
                        energys = regexp(thisfragmode,'energy" value="');
                        e = length(energys{1});
                        frag_value = cell(e,2);
                        for eng = 1:e
                            idx1 = energys{1}(eng)+15;
                            idx2 = refs{1}(eng)-1;
                            energy_text = thisfragmode{1}(idx1:idx2);
                            frag_value(eng,1) = {str2double(energy_text)};
                            ep = d + eng;
                            k1 = (names{1}(ep)+7);
                            k2 = (values{1}(ep))-1;
                            frag_type = thisfragmode{1}(k1:k2);
                            frag_value(eng,2) = {frag_type};
                        end
                        addinfosto{ii,jj} = frag_value;
                    catch
                        addinfosto{ii,jj} = []; % error('Error while reading file: Unable to convert ionization energy.');
                    end
                elseif strcmpi(addinfo{jj},'precision')
                    count = str2double(extractBetween(peaks_text,'count="','"'));
                    peaks_tag = regexp(peaks_text,'<binaryDataArray ');
                    peaks_close = regexp(peaks_text,'</binaryDataArray');
                    precision = cell(count,1);
                    for ll = 1:count
                        idx1 = peaks_tag{1}(ll);
                        idx2 = peaks_close{1}(ll);
                        binary_text = peaks_text{1}(idx1:idx2);
                        precision_text = extractBetween(binary_text,' name="','"');
                        precision = precision_text{1};
                    end
                    addinfosto{ii,jj} = precision;
                elseif strcmpi(addinfo{jj},'spectra')
                    if contains(peaks_text,'no compression')
                        count = str2double(extractBetween(peaks_text,'count="','"'));
                        peaks_tag = regexp(peaks_text,'<binaryDataArray ');
                        peaks_close = regexp(peaks_text,'</binaryDataArray');
                        try
                            binary = [];
                            for ll = 1:count
                                idx1 = peaks_tag{1}(ll);
                                idx2 = peaks_close{1}(ll);
                                binary_text = peaks_text{1}(idx1:idx2);
                                precision = extractBetween(binary_text,' name="','"');
                                binary_convert = extractBetween(binary_text,'<binary>','</binary>');
                                if strcmp(precision{1},'32-bit float')
                                    binary_uint = base64decode(binary_convert{1});
                                    binary_single = typecast(binary_uint,'single');
                                    binary(:,ll) = double(binary_single);
                                else
                                    binary_uint = base64decode(binary_convert{1});
                                    binary(:,ll) = typecast(binary_uint,'double');
                                end
                            end
                            addinfosto{ii,jj} = binary;
                        catch
                            addinfosto{ii,jj} = [];
                        end
                    else
                        error('Error: Cannot convert compressed binary data.');
                    end
                elseif strcmpi(addinfo{jj},'filterLine')
                    addinfosto{ii,jj} = extractBetween(scan_text,'"filter string" value="','"');
                end
            end
        end
        if mod(ii,200) == 0
            waitbar(ii/scancount,h,sprintf('%d/%d',ii,scancount));
        end
    end
    if ismember('std',infos)
        msdata.scannum = scanum;
        msdata.mslvl = mslvl;
        msdata.retime = retime;
        msdata.fragmode = fragmode;
        msdata.precursormz = precursormz;
        msdata.allprecursormz = allprecursormz;
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
end
close(h);
end
