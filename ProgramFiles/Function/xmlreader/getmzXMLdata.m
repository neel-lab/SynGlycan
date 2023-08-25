function msdata = getmzXMLdata(xmlobj,infos,h)
%GETMZXMLDATA: This function reads selected fields available in mzXML
% files.
%
% Syntax: 
%   msdata = getmzXMLdata(xmlobj,infos,h)
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
%   msdata = getmzXMLdata('a filename',{'std','spectra','totIonCurrent'},h = waitbar(0,'','Name','Writing msdata file...'));
%
%See also: GETMSDATA, GETINDEXEDMZMLDATA, GETMZMLDATA, PROTEOWIZARD
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
    manuinfoind = ismember(lower(manufacturer),lower(infos));
    if ~any(manuinfoind)
        manuinfo = lower(manufacturer{1});
    else
        manuinfo = manufacturer{manuinfoind};
    end  
    scan_tag = regexp(xmlobj,'scan num=');
    tagClose = regexp(xmlobj,'</scan');
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
    msdata.parentFile = extractBetween(xmlobj,'fileName="','fileType');
    for ii = 1:scancount
        scan_text = xmlobj((scan_tag(ii)):(tagClose(ii)-1));
        precursor_text = extractBetween(scan_text,'<precursorMz ','</');
        peaks_text = extractBetween(scan_text,'<peaks ','</');
        if ismember('std',infos)
            scanum(ii) = str2double(extractBetween(scan_text,'scan num="','"'));
            mslvl(ii) = str2double(extractBetween(scan_text,'msLevel="','"'));
            A = extractBetween(scan_text,'retentionTime="','"');
            B = regexp(A,'\d*','Match');
            try
                C = strcat((B{1}{1}),'.',(B{1}{2}));
            catch
                C = B{1};
            end
            retime(ii) = str2double(C);
            if ~isempty(precursor_text)
                precursormz(ii) = str2double(extractAfter(precursor_text,'">'));
                try
                    charge(ii) = str2double(extractBetween(precursor_text,'Charge="','"'));
                catch
                    charge(ii) = 0;
                end
                thisfragmode = extractBetween(precursor_text,'activationMethod="','"');
                fragmode{ii} = thisfragmode;
            else
                precursormz(ii) = -2;
                fragmode{ii} = '';
                charge(ii) = -1;
            end
        end
        if ~isempty(stdinfo)
            for jj = 1:length(stdinfo)
                if strcmpi(stdinfo{jj},'fragmode')
                    fragmode = extractBeween(precursor_text,'activationMethod="','"');
                    stdinfosto{ii,jj} = fragmode;
                elseif strcmpi(stdinfo{jj},'retentionTime')
                    stdinfosto{ii,jj} = (extractBetween(scan_text,'retentionTime="','"'));
                elseif strcmpi(stdinfo{jj},'num')
                    stdinfosto{ii,jj} = str2double(extractBetween(scan_text,'scan num="','"'));
                elseif strcmpi(stdinfo{jj},'msLevel')
                    stdinfosto{ii,jj} = str2double(extractBetween(scan_text,'msLevel="','"'));
                elseif strcmpi(stdinfo{jj},'precursorMz')
                    if ~isempty(precursor_text)
                        try
                            stdinfosto{ii,jj} = str2double(extractAfter(precursor_text,'">'));
                        catch
                            stdinfosto{ii,jj} = 0;
                        end
                    else
                        stdinfosto{ii,jj} = -2;
                    end
                elseif strcmpi(stdinfo{jj},'precursorCharge')
                    if ~isempty(precursor_text)
                        try
                            stdinfosto{ii,jj} = str2double(extractBetween(precursor_text,'Charge="','"'));
                        catch
                            stdinfosto{ii,jj} = 0;
                        end
                    else
                        stdinfosto{ii,jj} = -1;
                    end
                end
            end
        end        
        if ~isempty(addinfo)
            for jj = 1:length(addinfo)
                if strcmpi(addinfo{jj},'peaksCount')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'peaksCount="','"'));
                elseif strcmpi(addinfo{jj},'polarity')
                    addinfosto{ii,jj} = (extractBetween(scan_text,'polarity="','"'));
                elseif strcmpi(addinfo{jj},'scanType')
                    addinfosto{ii,jj} = (extractBetween(scan_text,'scanType="','"'));
                elseif strcmpi(addinfo{jj},'centroided')
                    addinfosto{ii,jj} = (extractBetween(scan_text,'centroided="','"'));
                elseif strcmpi(addinfo{jj},'deisotoped')
                    addinfosto{ii,jj} = (extractBetween(scan_text,'deisotoped="','"'));
                elseif strcmpi(addinfo{jj},'chargeDeconvoluted')
                    addinfosto{ii,jj} = (extractBetween(scan_text,'chargeDeconvoluted="','"'));
%                 elseif strcmpi(addinfo{j},'startMz')
%                     addinfosto{i,j} = ;
%                 elseif strcmpi(addinfo{j},'endMz')
%                     addinfosto{i,j} = ;
                elseif strcmpi(addinfo{jj},'lowMz')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'lowMz="','"'));
                elseif strcmpi(addinfo{jj},'highMz')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'highMz="','"'));
                elseif strcmpi(addinfo{jj},'basePeakMz')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'basePeakMz="','"'));
                elseif strcmpi(addinfo{jj},'basePeakIntensity')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'basePeakIntensity="','"'));
                elseif strcmpi(addinfo{jj},'totIonCurrent')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'totIonCurrent="','"'));
                elseif strcmpi(addinfo{jj},'precursorScanNum')
                    try
                        addinfosto{ii,jj} = str2double(extractBetween(scan_text,'precursorScanNum="','"'));
                    catch
                        addinfosto{ii,jj} = 0;
                    end
                elseif strcmpi(addinfo{jj},'ionisationEnergy')
                    addinfosto{ii,jj} = str2double(extractBetween(scan_text,'collisionEnergy="','"'));
                elseif strcmpi(addinfo{jj},'precision')
                    addinfosto{ii,jj} = str2double(extractBetween(peaks_text,'precision="','"'));
%                 elseif strcmpi(addinfo{j},'spectra')
%                     compression_text = extractBetween(scan_text,'compressionType="','"');
%                     if strcmpi(compression_text,'none')
%                         try
%                             binary = [];
%                                 binary_text = extractBetween(scan_text,'"m/z-int">','</peaks>');
%                                 precision = extractBetween(scan_text,'precision="','"');
%                                 if strcmp(precision{1},'32')
%                                     binary_uint = base64decode(binary_text{1});
%                                     binary_single = typecast(binary_uint,'single');
%                                     binary = double(binary_single);
%                                 else
%                                     binary_uint = base64decode(binary_text{1});
%                                     binary = typecast(binary_uint,'double');
%                                 end
%                             addinfosto{i,j} = binary;
%                         catch
%                             addinfosto{i,j} = []; 
%                         end
%                     else
%                         addinfosto{i,j} = [];
%                     end
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