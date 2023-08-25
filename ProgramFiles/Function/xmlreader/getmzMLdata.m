function msdata = getmzMLdata(xmlobj,infos,h)
%GETMZMLDATA: This function reads selected fields available in mzML files.
%
% Syntax:
% msdata = getmzMLdata(xmlobj,infos,h)
%
% Input:
% xmlobj: n x 1 cell array of strings. Expt. data file content in XML
%     format.
% infos: 1 x n cell array of strings. Name of the items to be extracted.
% h: handle. Waitbar.
%
% Output:
% msdata: structure. Extracted inforamtion. See GETMSDATA for detail.
%
% Examples:
% h = waitbar(0,'','Name','Writing msdata file...');
% msdata = getmzMLdata("a filename",{'std','spectra','totIonCurrent'},h);
%
% See also:
% GETMSDATA  GETINDEXEDMZMLDATA  GETMZXMLDATA  PROTEOWIZARD
%

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 11/10/2020

% Field names available. These names are similar to the names used in .mzML
%    format. If input "infors" contains 'std' all 6 fields will be extracted.
standardinfo = {'num','msLevel','retentionTime',...
    'precursorMz','precursorCharge','fragmode'};
additionalinfo = {'peaksCount','polarity',...
    'centroided','deisotoped','chargeDeconvoluted',...
    'startMz','endMz','lowMz',...
    'highMz','basePeakMz','basePeakIntensity',...
    'totIonCurrent','precursorScanNum','ionisationEnergy',...
    'precision','spectra','filterLine'}; %'scanType'
manufacturer = {'Thermo'};
if nargin == 0  % No input - Show available inputs as a form of help document
    disp(strjoin(['Standard: ',standardinfo],', '))
    disp(strjoin(['Custom: ',additionalinfo],', '))
    disp(['Default MS inst. manufacturer: ',manufacturer{1}])
    msdata = [];
else
    if isempty(infos)  % By default standard package is used
        infos = {'std'};
    elseif ischar(infos)
        % If input is string, only 1 field or standard package will be
        %     extracted
        infos = {infos};
    end
    
    % Fixing inputs uppercase/lowercase issue by converting everything to
    %     lowercase. Also only supported fields will be extracted
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
        charge = zeros(size(scanum));
    end
    if ~isempty(stdinfo)
        stdinfosto = cell(scancount,length(stdinfo));
    end
    if ~isempty(addinfo)
        addinfosto = cell(scancount,length(addinfo));
    end
    for ii = 1:scancount
        scan_text = xmlobj((scan_tag(ii)):(tagClose(ii)-1));
        if ~contains(scan_text,'"spectrum with no data"')
            scanlist_text = extractBetween(scan_text,'<scanList ','</scanList');
            precursor_text = extractBetween(scan_text,'<precursorList ','</precursorList');
            peaks_text = extractBetween(scan_text,'<binaryDataArrayList ','</binaryDataArrayList');
            if ismember('std',infos)
                scanum(ii) = str2double(extractBetween(scan_text,'id="scan=','"'));
                mslvl(ii) = str2double(extractBetween(scan_text,'ms level" value="','"'));
                A = extractBetween(scan_text,'scan time" value="','"');
                B = regexp(A,'\d*','Match');
                try  % Retention time - try to cover as many formats as possible
                    C = strcat((B{1}{1}),'.',(B{1}{2}));
                catch
                    C = B{1};
                end
                retime(ii) = str2double(C);
                if ~isempty(scanlist_text)
                    precursormz(ii) = str2double(extractBetween(scanlist_text,'Monoisotopic M/Z:" value="','"'));
                end
                if ~isempty(precursor_text)  % Tandem MS scans
                    try
                        charge(ii) = str2double(extractBetween(precursor_text,'charge state" value="','"'));
                    catch
                        charge(ii) = -1;
                    end
                    thisfragmode = extractBetween(precursor_text,'<activation>','</activation>');
                    names = regexp(thisfragmode,' name="');
                    dissociations = regexp(thisfragmode,' dissociation');
                    frag_count = length(dissociations{1});
                    this = thisfragmode{1};
                    if frag_count == 1  % Simple fragmentation method - CID HCD ETD
                        idx1 = names{1} + 7;  % Fixed format - directly extracting fragmentation method. 
                        idx2 = dissociations{1};
                        frag_text = this(idx1:idx2);
                        if strcmpi(frag_text,'collision-induced ')
                            fragmode{ii} = 'CID';
                        elseif strcmpi(frag_text,'beam-type collision-induced ')
                            fragmode{ii} = 'HCD';
                        elseif strcmpi(frag_text,'electron transfer ')
                            fragmode{ii} = 'ETD';
                        end
                    else  % Compounded fragmentation method - ETCID ETHCD
                        for ff = 1:frag_count
                            idx1 = (names{1}(ff))+7;
                            idx2 = dissociations{1}(ff);
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
                else  % MS1 scans
                    precursormz(ii) = -2;
                    fragmode{ii} = '';
                    charge(ii) = -1;
                end
            end
            if ~isempty(stdinfo)  % Extract information defined in "standardinfo"
                for jj = 1:length(stdinfo)
                    if strcmpi(stdinfo{jj},'fragmode')
                        thisfragmode = extractBetween(precursor_text,'<activation>','</activation>');
                        names = regexp(thisfragmode,' name="');
                        dissociations = regexp(thisfragmode,' dissociation');
                        frag_count = length(dissociations{1});
                        this = thisfragmode{1};
                        if frag_count == 1  % Simple fragmentation
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
                        else  % Compounded fragmentation
                            for ff = 1:frag_count
                                idx1 = (names{1}(ff))+7;
                                idx2 = dissociations{1}(ff);
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
                        stdinfosto{ii,jj} = (extractBetween(scan_text,'retentionTime="','"'));
                    elseif strcmpi(stdinfo{jj},'num')
                        stdinfosto{ii,jj} = str2double(extractBetween(scan_text,'scan num="','"'));
                    elseif strcmpi(stdinfo{jj},'msLevel')
                        stdinfosto{ii,jj} = str2double(str2double(extractBetween(scan_text,'msLevel="','"')));
                    elseif strcmpi(stdinfo{jj},'precursorMz')
                        if ~isempty(precursor_text)
                            try
                                stdinfosto{ii,jj} = str2double(extractBetween(precursor_text,'"m/z" value="','"'));
                            catch
                                stdinfosto{ii,jj} = 0;
                            end
                        else
                            stdinfosto{ii,jj} = -2;
                        end
                    elseif strcmpi(stdinfo{jj},'precursorCharge')
                        if ~isempty(precursor_text)
                            try
                                stdinfosto{ii,jj} = str2double(extractBetween(precursor_text,'charge state" value="','"'));
                            catch
                                stdinfosto{ii,jj} = 0;
                            end
                        else
                            stdinfosto{ii,jj} = -1;
                        end
                    end
                end
            end
            if ~isempty(addinfo)  % Extract information defined in "additionalinfo"
                for jj = 1:length(addinfo)
                    addinfosto{ii,jj} = -1;
                    if strcmpi(addinfo{jj},'peaksCount')
                        addinfosto{ii,jj} = str2double(extractBetween(scan_text,'defaultArrayLength="','"'));
                    elseif strcmpi(addinfo{jj},'polarity')
                        if contains(xmlobj,'positive scan')
                            addinfosto{ii,jj} = '+';
                        elseif contains(xmlobj,'negative scan')
                            addinfosto{ii,jj} = '-';
                        else
                            addinfosto{ii,jj} = [];
                        end
                    elseif strcmpi(addinfo{jj},'centroided')
                        if contains(xmlobj,'centroid mass spectrum')
                            addinfosto{ii,jj} = 'centroided';
                        else
                            addinfosto{ii,jj} = [];
                        end
                    elseif strcmpi(addinfo{jj},'deisotoped')
                        addinfosto{ii,jj} = extractBetween(xmlobj,'"deisotoping" value="','"');
                    elseif strcmpi(addinfo{jj},'chargeDeconvoluted')
                        addinfosto{ii,jj} = extractBetween(xmlobj,'"charge deconvolution" value="','"');
                    elseif strcmpi(addinfo{jj},'startMz')
                        addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"scan window lower limit" value="','"'));
                    elseif strcmpi(addinfo{jj},'endMz')
                        addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"scan window upper limit" value="','"'));
                    elseif strcmpi(addinfo{jj},'lowMz')
                        addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"lowest m/z value" value="','"'));
                    elseif strcmpi(addinfo{jj},'highMz')
                        addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"highest m/z value" value="','"'));
                    elseif strcmpi(addinfo{jj},'basePeakMz')
                        addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"base peak m/z" value="','"'));
                    elseif strcmpi(addinfo{jj},'basePeakIntensity')
                        addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"base peak intensity" value="','"'));
                    elseif strcmpi(addinfo{jj},'totIonCurrent')
                        addinfosto{ii,jj} = str2double(extractBetween(scan_text,'"total ion current" value="','"'));
                    elseif strcmpi(addinfo{jj},'precursorScanNum')
                        if any(strfind(scan_text,'master scan number'))
                            temp = str2double(extractBetween(scan_text,'name="master scan number" value="','"'));
                            if ~isempty(temp)
                                addinfosto{ii,jj} = temp;
                            end
                        else
                            try
                                temp = str2double(extractBetween(precursor_text,'"scan=','"'));
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
                            if isempty(thisfragmode)
                                addinfosto{ii,jj} = '';
                            else
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
                            end
                        catch
                            error('Error while reading file: Unable to convert ionization energy.'); % addinfosto{i,j} = [];
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
                            error('Cannot convert compressed binary data.');
                        end
                    elseif strcmpi(addinfo{jj},'filterLine')
                        addinfosto{ii,jj} = extractBetween(scan_text,'"filter string" value="','"');
                    end
                end
            end
        else
            scanum(ii) = str2double(extractBetween(scan_text,'id="scan=','"'));
            mslvl(ii) = str2double(extractBetween(scan_text,'ms level" value="','"'));
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