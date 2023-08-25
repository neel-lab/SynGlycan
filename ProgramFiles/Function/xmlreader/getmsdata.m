function msdata = getmsdata(filename,infos,intermediatefileloc)
%GETMSDATA: This function is the front end for converting a raw or XML file
% containing MS data into a .mat file to use for analysis
%
% Syntax:
% msdata = getmsdata(filename,infos,intermediatefileloc)
%
% Input:
% filename: string. Full path to experiment data file.
% infos: 1 x n cell array of strings. Information to be read and saved.
%     Available items (the strings in cell) are:
%     'std': standard information package, usually stored in scan header. These
%         includes:
%             'num': Scan number
%             'msLevel': MS level
%             'retentionTime': Retention time
%             'precursorMz': Precursor ion's m/z
%             'precursorCharge': Precursor ion's charge
%             'fragmode': Fragmentation method
%     'peaksCount': number of peaks in this spectrum
%     'polarity': positive or negative scan
%     'scanType': scan type
%     'centroided': peaks are centroided or not (profile mode)
%     'deisotoped': peaks are de-isotoped or not
%     'chargeDeconvoluted': peaks have their chage deconvoluted or not
%     'startMz': scan window lower limit (m/z)
%     'endMz': scan window upper limit (m/z)
%     'lowMz': min. m/z in spectrum
%     'highMz': max. m/z in spectrum
%     'basePeakMz': m/z of base peak
%     'basePeakIntensity': intensity of base peak
%     'totIonCurrent': total ion current
%     'precursorScanNum': scan number of precursor scan
%     'ionisationEnergy': ionisation energy
%     'precision': precision
%     'spectra': spectrum m/z - intensity list
%     'filterLine':  (Thermo specific) information in section "filterLine"
% intermediatefileloc: string. The location to store intermediate files.
%
% Output:
% msdata: structure. Information extracted from file.
%
% Note:
% Currently only .mzML , .RAW and .mat format is supported. Support
%     for other formats will be available in future release.
%
% Examples:
% Users are encouraged to test this function with their own data.
% It is recommended to use raw data files that can be converted to .mzML
%     files by Proteowizard.
% For input "infos" this value is recommended:
%     {'std','totIonCurrent','precursorScanNum','spectra'}
%
% See also:
% PREPROCESSGUI  GETINDEXEDMZMLDATA  GETMZMLDATA  GETMZXMLDATA  PROTEOWIZARD
%

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 11/10/2020

[msfilepath,name,ext] = fileparts(filename);  % Find file type by its extension
if strcmpi(ext,'.RAW')
    if ~strcmpi(intermediatefileloc(end),'\')
        intermediatefileloc = [msfilepath,'\'];
    end
    [filepath,~,~] = fileparts(intermediatefileloc);
    [status,cmdout] = proteowizard(filename,filepath); %switched from filename to name
    if status == 0
        filename = [filepath,'\',name,'.mzML'];
    else
        error('Error while converting .raw file.')
        disp(cmdout);
    end
end
j = msgbox('Reading XML file...');
xmlobj = fileread(filename);  % Read the entire XML file
close(j);
h = waitbar(0,'','Name','Writing msdata file...');
breaks = regexp(xmlobj, '[\n]');
first_line = xmlobj(1:(breaks(1)-1));  % Read file header section

% Function "EXTRACTBETWEEN" costs memory space, here instead of putting in
%     the entire file, we put in a large enough section only
spaces = extractBetween(xmlobj(1:breaks(min(20,length(breaks)))-1),'<',' ');
file_type = spaces{2};
if strcmpi(file_type,'mzML')
    msdata = getmzMLdata(xmlobj,infos,h);
elseif strcmpi(file_type,'indexedmzML')
    msdata = getIndexedmzMLdata(xmlobj,infos,h);
elseif strcmpi(file_type,'mzXML')
    msdata = getmzXMLdata(xmlobj,infos,h);
elseif strcmpi(first_line,'BEGIN IONS')
    msdata = XgetMGFdata(xmlobj,infos,h);
end
end
