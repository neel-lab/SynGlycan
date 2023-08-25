function scoreCSVwrite(outputfilefullname,handles)
% SCORECSVWRITE: Write the scoring results to a csv file
%
% Syntax:
% scoreCSVwrite(csvheader,datastructure,outputfilefullname)
%
% Input:
% csvheader: structure. Header of result CSV file. These information will
%     be written at the beginning of the file. This function has already
%     defined the format of certain fields and will ignore other
%     informations.
% datastructure: structure. Scoring results. Each field in
%     "datastructure" will occupy a column. The sequence of these columns
%     remains the same.
% outputfilefullname: String. Full file name for the output file.
%
% Output:
% A .csv file containing header and scoring results.
%
% Note:
% Field names in "csvheader" must be all lowercase.
% 
% Example:
% N/A
% 
% Children function:
% N/A
% 
% See Also:
% N/A

% Author: Gang Liu, rewritten by Kai Cheng
% Date Lastly Updated: 01/01/2020

% write header to output file
scoredata = handles.scoredata;
displayoptions = handles.displayoptions;
scoreintdata = displayoptions.scoreintdata;
scoreoptions = displayoptions.scoreintdata.scoreoptions;
csvheader.PeptideFile = displayoptions.scoreintdata.sliminput.pepfile;
csvheader.ExperimentData = displayoptions.scoreintdata.sliminput.exptdata;
switch scoreoptions.isHCDtrigger
    case true
        csvheader.IsHCDTrigger = 'Yes';
    case false
        csvheader.IsHCDTrigger = 'No';
end
csvheader.MS1Tolerence = [num2str(scoreoptions.ms1tol),' ',...
    scoreoptions.ms1tolunit];
MS2Tolerence = '';
for ii = 1:length(scoreoptions.analyzefragmode)
    MS2Tolerence = [MS2Tolerence,scoreoptions.analyzefragmode{ii},': ',...
        num2str(scoreoptions.ms2tol(ii)),' ',...
    scoreoptions.ms2tolunit{ii},'; '];
end
csvheader.MS2Tolerence = MS2Tolerence;
FragmentationNumber = 'NpFrag/NgFrag/NpFrag: ';
for ii = 1:length(scoreoptions.analyzefragmode)
    FragmentationNumber = [FragmentationNumber,'(',scoreoptions.analyzefragmode{ii},') ',...
        num2str(scoreoptions.fragnum(ii,:)),'; '];
end
csvheader.FragmentationNumber = FragmentationNumber;
LabileMonosaccharide = 'Neu5Ac and Fuc are labile: ';
for ii = 1:length(scoreoptions.analyzefragmode)
    if scoreoptions.monosacislabile(ii)
        LabileMonosaccharide = [LabileMonosaccharide,'(',scoreoptions.analyzefragmode{ii},') ',...
            'Yes; '];
    else
        LabileMonosaccharide = [LabileMonosaccharide,'(',scoreoptions.analyzefragmode{ii},') ',...
            'No; '];
    end
end
csvheader.LabileMonosaccharide = LabileMonosaccharide;
GlycanPeptideSimultaneousFragmentation = 'Glycan and peptide may fragment simultaneously: ';
for ii = 1:length(scoreoptions.analyzefragmode)
    if scoreoptions.simultaneousfrag(ii)
        GlycanPeptideSimultaneousFragmentation = [GlycanPeptideSimultaneousFragmentation,...
            '(',scoreoptions.analyzefragmode{ii},') ','Yes; '];
    else
        GlycanPeptideSimultaneousFragmentation = [GlycanPeptideSimultaneousFragmentation,...
            '(',scoreoptions.analyzefragmode{ii},') ','No; '];
    end
end
csvheader.GlycanPeptideSimultaneousFragmentation = GlycanPeptideSimultaneousFragmentation;
AdditionalOxoniumIon = 'Add additional oxonium fragment ions';
for ii = 1:length(scoreoptions.analyzefragmode)
    if scoreoptions.addoxoniumion(ii)
        AdditionalOxoniumIon = [AdditionalOxoniumIon,...
            '(',scoreoptions.analyzefragmode{ii},') ','Yes; '];
    else
        AdditionalOxoniumIon = [AdditionalOxoniumIon,...
            '(',scoreoptions.analyzefragmode{ii},') ','No; '];
    end
end
csvheader.AdditionalOxoniumIon = AdditionalOxoniumIon;
PeptideFragmentType = '';
for ii = 1:length(scoreoptions.analyzefragmode)
    PeptideFragmentType = [PeptideFragmentType,scoreoptions.analyzefragmode{ii},': ',...
        scoreoptions.userpepiontyp{ii},'; '];
end
csvheader.PeptideFragmentType = PeptideFragmentType;
csvheader.MaximumGlycanStubLength = num2str(scoreoptions.maxstublen);
csvheader.CutOffMedian = num2str(scoreoptions.cutoffmed);
csvheader.FracMax = num2str(scoreoptions.fracmax);
if isempty(scoreoptions.selectpeak)
    csvheader.SelectivePeak = 'None';
else
    csvheader.SelectivePeak = num2str(scoreoptions.selectpeak);
end
displaydataind = handles.displaydataind;
usercustom = handles.usercustom;
if iscell(displaydataind)
    writedisplaydataind = [];
    for i = 1:size(displaydataind,1)
        writedisplaydataind = [writedisplaydataind;cell2mat(displaydataind(i,:))];
    end
elseif isnumeric(displaydataind)
    writedisplaydataind = displaydataind;    
end
write_displayoptions = displayoptions;
write_displayoptions.combineisomer = false;
[writedata,~,displaycolnames,~] = ...
    builddisptable(scoredata,writedisplaydataind,write_displayoptions,usercustom);
displaycolnames = strrep(displaycolnames,'%','Percentage_');
for ii = 1:length(displaycolnames)
    outputscoredata.(displaycolnames{ii}) = writedata(:,ii);
end

% fldnms = fieldnames(csvheader);
% for ii = 1:length(fldnms)
%     thisfldnm = fldnms{ii};
%     fprintf(FID,'%s\n',thisfldnm);
%     fprintf(FID,'%s\n',csvheader.(thisfldnm));
% end
% struct2csv(datastructure,FID);  % Write entire structure to file including one line headerend
% fclose(FID);

datafldnms = fieldnames(outputscoredata);
tblheight = length(outputscoredata.(datafldnms{1}));
firstsheetdata2 = cell(tblheight + 1,length(datafldnms));
firstsheetdata2(1,:) = datafldnms;
for ii = 1:length(datafldnms)
    firstsheetdata2(2:end,ii) = outputscoredata.(datafldnms{ii});
end
firstsheetdata = cell(length(csvheader),length(datafldnms));
headerfldnms = fieldnames(csvheader);
for ii = 1:length(headerfldnms)
    firstsheetdata{ii*2 - 1,1} = headerfldnms{ii};
    firstsheetdata{ii*2,1} = csvheader.(headerfldnms{ii});
end
firstsheetdata = [firstsheetdata;firstsheetdata2];
blocksize = 5000;
numblocks = ceil(size(firstsheetdata,1)/blocksize);
for ii = 1:numblocks
    writecell(firstsheetdata(((ii - 1)*blocksize + 1):min(size(firstsheetdata,1),ii*blocksize),:),outputfilefullname,'WriteMode','append','Sheet','GlycopeptidePSMs');
end

if scoreintdata.scoreoptions.isHCDtrigger
    summarydata = scoredata(writedisplaydataind(:,1));
    proteins = {summarydata.Protein};
    sgps = {summarydata.SGP};
    peptides = cell(size(summarydata));
    glycans = cell(size(summarydata));
    glypos = cell(size(summarydata));
    scans = [summarydata.Scan];
    quants = [summarydata.Quant];
    pepstartpos = zeros(size(summarydata));
    [~,~,ind] = unique(scans);
    for ii = 1:max(ind)
        tempquant = quants(ind == ii);
        tempquant = tempquant/length(tempquant);
        quants(ind == ii) = tempquant;
    end
    for ii = 1:length(sgps)
        [p,g,~] = breakGlyPep(sgps{ii});
        peptides{ii} = p.pep;
        glycans{ii} = {g.struct};
        glypos{ii} = [g.pos];
    end
    if exist(scoreintdata.sliminput.pepfile,'file')
        pepfileloc = scoreintdata.sliminput.pepfile;
    else
        [f,p] = uigetfile('*.txt','Locate the digestion file');
        pepfileloc = fullfile(p,f);
    end
    digestfileinfo = digestfileanaly(pepfileloc,1);
    
    if isempty(digestfileinfo.protseqpath)
        [f,p] = uigetfile('*.txt','Locate the FASTA file');
        fastafileloc = fullfile(p,f);
    else
        fastafileloc = digestfileinfo.protseqpath{1};
    end
    if ~exist(fastafileloc,'file')
        [f,p] = uigetfile('*.txt','Locate the FASTA file');
        fastafileloc = fullfile(p,f);
    end
    protseq = fastaread(fastafileloc);
    fastaprotheaders={protseq.header};
    summaryprotheaders={summarydata.Protein};
    fastaseq='';
    for ii = 1:length(summarydata) %for each glycan's peptide
        protstr=summaryprotheaders(ii); %where is this protein located in the fasta, character array
        for jj = 1:length(fastaprotheaders)
            fastastr=fastaprotheaders(jj);
            comparison=strcmpi(protstr,fastastr);
            if comparison == 1 %If protein from summary is protein in fasta
                fastaseq = protseq(jj).sequence; %get corresponding prot seq
%                 pepstartpos(ii) = strfind(fastaseq,peptides{ii}); %find the start position of this peptide sequence in this protein and save it
            end
        end
    end
    
    %% SHEET: PROTEIN
    data_proteinsheet = {'Prot. Name','Prot. #PSM','Prot. AUC','Peptides',...  % 1 2 3 4
        'Pep. #PSM','Gly. Pos.','Glycans','Pep. AUC'};  % 5 6 7 8
    [uniprot,~,uniprotind] = unique(proteins);
    for ii = 1:max(uniprotind)
        thisprot_proteinsheet = {};
        thisuniprot = uniprot{ii};
        tempsummarydata = summarydata(uniprotind == ii);
        uniprot_peptides = peptides(uniprotind == ii);
        uniprot_glycans = glycans(uniprotind == ii);
        uniprot_glycopos = glypos(uniprotind == ii);
        uniprot_pepstartpos = pepstartpos(uniprotind == ii);
        uniprot_quant = quants(uniprotind == ii);
        [unipep,~,unipepind] = unique(uniprot_peptides);
        for jj = 1:max(unipepind)
            tempthispep_uniprot_glycopos = uniprot_glycopos(unipepind == jj);
            thispep_uniprot_glycopos = [];
            tempthispep_uniprot_glycans = uniprot_glycans(unipepind == jj);
            thispep_uniprot_glycans = {};
            for kk = 1:length(tempthispep_uniprot_glycopos)
                thispep_uniprot_glycopos = [thispep_uniprot_glycopos,tempthispep_uniprot_glycopos{kk}];
                thispep_uniprot_glycans = [thispep_uniprot_glycans,tempthispep_uniprot_glycans{kk}];
            end
            thispep_uniprot_glycopos = thispep_uniprot_glycopos + uniprot_pepstartpos(1);
            thispep_uniprot_glycans = unique(thispep_uniprot_glycans);
            thispepdata_proteinsheet = cell(length(thispep_uniprot_glycans),length(data_proteinsheet));
            thispepdata_proteinsheet{1,4} = unipep{jj};
            thispepdata_proteinsheet{1,5} = sum(unipepind == jj);
            thispepdata_proteinsheet{1,6} = num2str(unique(thispep_uniprot_glycopos));
            thispepdata_proteinsheet(:,7) = thispep_uniprot_glycans;
            thispepdata_proteinsheet{1,8} = sum(uniprot_quant(unipepind == jj));
            thisprot_proteinsheet = [thisprot_proteinsheet;thispepdata_proteinsheet];
        end
        thisprot_proteinsheet(1,1:3) = {thisuniprot,length(tempsummarydata),...
            sum(uniprot_quant)};
        data_proteinsheet = [data_proteinsheet;thisprot_proteinsheet];
    end
    writecell(data_proteinsheet,outputfilefullname,'Sheet','Protein');
    
    %% SHEET: PEPTIDE
    
    
    %% SHEET: GLYCAN
end

























