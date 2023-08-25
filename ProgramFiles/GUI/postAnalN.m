function postAnalN(pepList,glyList,pData,gpData,OutputFolder,Out_N_Loc)
% Inputs:
% pepfile: Contains list of peptides in search [all peptides contain NXS/T sites, with some repetition due to missed cleavage]
% Nglycanfile: List of Nglycans in search
% gpDataFile: glycopeptide best isomer file
% pDataFile: 1_1 file for peptide HCD search

% Outputs:
% AUC2: AUC for each glycopepide and peptide structure
% Count2: Count for each glycopeptide and peptide structure
% AUC_Ucomp: AUC for each glycan composition
% Count_Ucomp: Count for each glycan composition
% AUC_classify: AUC for each glycan classification
% Count_classify: Count for each glycan classification

AUC=zeros(length(pepList),1+size(glyList,1));     % to store AUC for each peptide/glycopeptide
Scan=cell(length(pepList),1+size(glyList,1));     % to store AUC for each peptide/glycopeptide
Count=zeros(length(pepList),1+size(glyList,1));   % to count number of spectra for each (glyco)peptide
[gpDataP,gpDataG,gpDataMod]=cellfun(@(x) breakGlyPep(x), gpData.SGP, 'UniformOutput', false);
gpDataP=cellfun(@(x) x.pep, gpDataP, 'UniformOutput', false);  % extracts the peptide from the structure
gpDataG=cellfun(@(x) x.struct, gpDataG, 'UniformOutput', false);  % extracts the glycan structure from the structure
pQuant = pData.Quant;
gpQuant = gpData.AUC;

% Find the average AUC for a given HCD scan in gp data and equally divide
% among hits from given HCD scan
gpHCD = gpData.Scan_HCD;
[~, ~, c] = unique(gpHCD);
for i = 1:max(c)
    idx = find(c==i);
    avg=mean(gpQuant(idx))/length(idx);
    gpQuant(idx)=avg;
end

for i = 1:size(pepList,1)

    % first get the AUC for the peptide without glycan and sum that
    idx = cellfun(@(x) strcmp(pepList{i}, x), pData.SGP, 'UniformOutput', false);
    pQuantTemp=pQuant;
    pQuantTemp(~cell2mat(idx)) = 0;
    header{1} = '';
    AUC(i,1) = sum(pQuantTemp); %For pep data, sum AUC for given 1/12 peptide found in pep file and place in first column
    Scan(i,1) = {pData.Scan(find(cell2mat(idx)>0))}; %For pep data, put all scans containing given peptide in first column
    Count(i,1) = length(find(cell2mat(idx))); % For pep data, count how many 1/12 peptide hits there are and place in first column

    % second count the number of hits/spectra contributing to peptide AUC
    Count(i,1) = sum(cell2mat(idx)>0);
    parfor j = 1:size(glyList,1)   % parfor if you want parallel computing, else just use for

        % Third, get the AUC for individual glycopeptides and sum them
        header{1+j} = glyList{j};
        idx2 = cellfun(@(x) strcmp(pepList{i}, x), gpDataP, 'UniformOutput', false);              % Find all given instances of this peptide in gpData
        idx3 = cellfun(@(x) strcmp(glyList{j}, x), gpDataG, 'UniformOutput', false);              % Find all given instances of this glycan in gpData
        idx_combined = or(~cell2mat(idx2),~cell2mat(idx3));                                       
        gpQuantTemp = gpQuant;
        gpQuantTemp(idx_combined) = 0;                                                            % Remaining AUC in list goes to 0
        AUC(i,1+j) = sum(gpQuantTemp);                                                            
        Scan(i,1+j) = {gpData.Scan_HCD(idx_combined==0)};
        % Fourth, get number of spectra for each glycopeptide
        Count(i,1+j) = sum(idx_combined==0);

        % Tables will look like... (1007x12)

        % AUC table
        % Header      N-glycan1             N-glycan2 ...
        % Peptide1    AUC for g1p1          AUC for g2p1
        % Peptide2    AUC for g1p2          AUC for g2p2
        % ...         ...                   ...
        % ...         ...                   ...

        % Scan table
        % Header      N-glycan1                          N-glycan2 ...
        % Peptide1    HCD scans containing g1p1          HCD scans containing g2p1
        % Peptide2    HCD scans containing g1p2          HCD scans containing g2p2
        % ...         ...                                ...
        % ...         ...                                ...

        % Count table
        % Header      N-glycan1                               N-glycan2 ...
        % Peptide1    # of HCD scans containing g1p1          # of HCD scans containing g2p1
        % Peptide2    # of HCD scans containing g1p2          # of HCD scans containing g2p2
        % ...         ...                                     ...
        % ...         ...                                     ...
    end
end

% Get unique glycan composition, list of 1007 compositions including blank
% at first position
mono={'h','n','f','s'};
comp = cell(length(header),1);
for j = 1: length(header)
    temp = [];
    for i = 1:length(mono)
        temp = [temp,mono{i},num2str(count(header{j},mono{i}))];
    end
    comp{j}=temp;
end
Ucomp = unique(comp);

% build glycan classifier

% 1. high_mannose; 2. pauci_mannose; 3. hybrid; 4. complex; 5. fucosylated
% 6. core_fucosylated_only; 7. terminal_fucosylated_only; 8. both_core_and_terminal_fucosylated
% 9. terminated_by_Neu5Ac;  10. terminated_by_Gal; 11. terminated_by_GlcNAc;
% 12. mono_antennary; 13. bi_antennary; 14. tri_antennary = 0; 15. tetra_antennary;
% 16. other_antennary; 17.bisecting; 18. Lewis; 19. sialylLewis; 20. LacdiNAc
nClassifier = 20;                   % This is only done for glycopeptides (not peptides alone)
classifier = zeros(nClassifier,length(header)); %20 x 1007
for j = 2:length(header)
    classifier(:,j) = getglycanstructfeatures_2(header{j})';
    % classifier table
    % Header      N-glycan1    N-glycan2 ...
    % Class1      Y/N          Y/N
    % Class2      Y/N          Y/N
    % ...         ...          ...
    % Class20     ...   
end

% Get to unique peptide by removing items with missed cleavage
pep={};
for j = 1:length(pepList) %12 peptides
    str = string(pepList{j}); %add string
    arra = regexp(str,'[KR][^P]');
    if isempty(arra)
        pep=[pep,str];
    end
end


% All rows with identical pep are collapsed. This results in AUC2 and Count2
% This section also calculates AUC and Count matrices based on composition and glycan
% classifier
AUC2 = [];
Count2= [];
Scan2 = {};
AUC_Ucomp = zeros(length(pep),length(Ucomp)); %4 peps x 242 unique Ngly compositions
Count_Ucomp = zeros(length(pep),length(Ucomp)); %4 peps x 242 unique Ngly compositions
AUC_classify = zeros(length(pep),size(classifier,1)); %4 peps x 20 classes
Count_classify = zeros(length(pep),size(classifier,1)); %4 peps x 20 classes
for j = 1:length(pep)       % compress and remove duplicate peptides (done for peptides and glycopeptides)
    idx2 = cellfun(@(x) contains(x,pep{j}), pepList, 'UniformOutput', false); %1 peptide + 2 miscleaved peptides => 1 peptides only, repeat 4x
    idx2 = (cell2mat(idx2)>0);
    AUC2 = [AUC2;sum(AUC(idx2,:))]; %Get total AUC for glycans on given 1/4 peptides
    Count2 = [Count2;sum(Count(idx2,:))]; %Get total number of HCD scans for glycans on given peptide
    temp = cell(1,size(Scan,2)); %1x1007 glycans
    idx2 = find(idx2);
    for t = 1:length(idx2)  %1x1007, 12 x 1007
        temp = cellfun(@(x, y) [x;y], temp,Scan(idx2(t),:),'UniformOutput',false); %combines HCD scans from 3 peps x 1007 glycans into 1 pep x 1007 glycans
    end
    Scan2 = [Scan2;temp]; % 1 pep (no miscleavages) x 1007 glycans, HCD scans grouped per glycan
            % Scan2 table
        % Header      N-glycan1                          N-glycan2 ...
        % Peptide1    HCD scans containing g1p1          HCD scans containing g2p1
        % Peptide2    HCD scans containing g1p2          HCD scans containing g2p2
        % Peptide3    HCD scans containing g1p3          ...
        % Peptide4    HCD scans containing g1p4          ...

    for k = 1:length(Ucomp) % Ucomp included peptide case since 'n0h0s0f0' is incuded
        idx3 = cellfun(@(x) strcmp(x,Ucomp{k}), comp, 'UniformOutput', false); %find location of unique composition in composition (1007)
        tempAUC2 = AUC2(j,:); %Get total AUC for glycans on given peptide
        tempAUC2 (cell2mat(idx3)<1) = 0;
        AUC_Ucomp(j,k) = sum(tempAUC2); %Put total AUC in 4 pep x 242 glycan composition table
        tempCount2 = Count2(j,:);
        tempCount2 (cell2mat(idx3)<1) = 0;
        Count_Ucomp(j,k) = sum(tempCount2); %Put number of HCD scans in 4 pep x 242 glycan composition table
    end
    for m =1:size(classifier,1)  % should be nCalssifier of them; only for glycopeptides (not peptides)
        AUC_classify(j,m) = sum(classifier(m,:).*AUC2(j,:)); %for given class multiply (1 glycan class x 1007 glycans) cell by AUC in (1 peptide x 1007 glycans) cell, put in (4 peps x 20 glycan classes) table
        Count_classify(j,m) = sum(classifier(m,:).*Count2(j,:)); %for given class multiply (1 glycan class x 1007 glycans) cell by Counts of HCD hits in (1 peptide x 1007 glycans) cell, put in (4 peps x 20 glycan classes) table
    end
end
AUC_classify=[AUC2(:,1), AUC_classify];             % add data for peptide now, None is methionine oxidation or pep alone? %4 peptides x (1 Total peptide AUC from pep Quant + AUC of 20 glycan classes)
Count_classify=[Count2(:,1), Count_classify];       % 4 peptides x (1 Total peptide Counts from pep Quant + Counts of 20 glycan classes)
AUC_classify=[AUC_classify;sum(AUC_classify,1)];    % last row is sum of column data
Count_classify=[Count_classify;sum(Count_classify,1)]; % last row is sum of column data

pep=[pep,'all']';
% PSM grouping
PSMs=sum(Count2,2); %4 peptides x 1 total counts
PSMsTotal=sum(PSMs,1); %Total counts (HCD scans)
PSMs=num2cell([PSMs;PSMsTotal]);

% glycopeptide grouping
temp_Count2=double(logical(Count2)); %4 peptides x 1007 glycans indicating which combinations contain hits
GP_1pep=sum(temp_Count2,2); % 4 peptides x 1 column, How many unique N-glycan structures were found on given peptide out of 1007
GP=num2cell([GP_1pep;sum(GP_1pep)]); %How many unique N-glycopeptide structures were found on given peptide, 5th row is total unique peptide/glycan combinations
%temp_GP=logical(sum(temp_Count2,1));
%GP_allpep=sum(double(temp_GP));
%GP=num2cell([GP_1pep;GP_allpep]);

% Glycan structure grouping
GlycanStruct_Array=double(logical(Count2)); %4 peptides x 1007 glycans indicating which combinations contain hits
GlycanStruct_1pep=sum(GlycanStruct_Array,2); %4 peptides x 1 column, How many unique N-glycan structures were found on given peptide out of 1007, # of unique glycopeptides
temp_GlycanStruct=logical(sum(GlycanStruct_Array,1)); %1 x 1007 glycans indicating which glycans contain hits
GlycanStruct_allpep=sum(double(temp_GlycanStruct)); % How many unique glycan structures out of 1007 contain hits
GlycanStruct=num2cell([GlycanStruct_1pep;GlycanStruct_allpep]); %4 peptides x 1 column, How many unique N-glycan structures were found on given peptide out of 1007, 
            % 5th row is how many unique glycan structures out of 1007 unique glycan structures were found

% Composition grouping
GlycanComp_Array=double(logical(Count_Ucomp)); %Which 4 pep x 242 unique glycan compositions are present, T/F
GlycanComp_1pep=sum(GlycanComp_Array,2); %4 peps x 1 column, How many unique compositions are present on given 1/4 peptide
temp_GlycanComp=logical(sum(GlycanComp_Array,1)); %1 row x 242 unique glycan compositions, Which compositions were found in whole sample
GlycanComp_allpep=sum(double(temp_GlycanComp)); % 1 x 1, How many unique compositions were found
GlycanComp=num2cell([GlycanComp_1pep;GlycanComp_allpep]); %5th row is how many unique glycan compositons out of 242 unique glycan compositions were found

% Normlize AUC classify
classes={'None','high_mannose','pauci','hybrid','complex','non-Fuc',...
    'coreFuc','termFuc','bothFuc','Neu5AcTerm','GalTerm','GlcNAcTerm','mono','bi','tri','tetra',...
    'other','bisect','non-bisect','lacdinac','lewis','slewis'};
Norm_AUC_classify = zeros(length(pep),size(classifier,1)); %(4 peptides + all) x 20 classes
for j = 1: size(Norm_AUC_classify,1)
    Norm_AUC_classify(j,1) = 100*AUC_classify(j,1)/sum(AUC_classify(j,1:21));       % all types of peptides
    Norm_AUC_classify(j,2:5) = 100*AUC_classify(j,2:5)/sum(AUC_classify(j,2:5));    % for all types of glycans
    Norm_AUC_classify(j,6:9) = 100*AUC_classify(j,6:9)/sum(AUC_classify(j,4:5));    % for all Fuc type (based on hybrid and complex)
    Norm_AUC_classify(j,10:18) = 100*AUC_classify(j,10:18)/AUC_classify(j,5);       % all others are based on complex
    Norm_AUC_classify(j,19) = 100-Norm_AUC_classify(j,18);                          % adding non-bisecting
    Norm_AUC_classify(j,20:22) = 100*AUC_classify(j,19:21)/AUC_classify(j,5);       % all others are based on complex
end

[infilepath,infilename,infileext] = fileparts(Out_N_Loc);
OutputFolder=char(OutputFolder);
outfilename = strcat(OutputFolder,'\0_',infilename,'.xlsx');
outfilename = strrep(outfilename,'.xlsx','_NGlyPostAnal.xlsx');

outputAUC=[cellstr(pep), PSMs, GP, GlycanStruct, GlycanComp,num2cell(Norm_AUC_classify)]; %(5 rows, peptides + all), # of HCD scans, # unique glycopeptide structures, # unique glycan structures out of 1007, # unique glycan compositions out of 242, Norm AUC
%outputAUC=[outputAUC(3:4,:);outputAUC(1:2,:);outputAUC(5,:)];
outputAUC_Table=cell2table(outputAUC);
h1=["Peptide", "PSM", "Unique GP","Unique GlyStruct", "Unique Comp", classes];
outputAUC_Table.Properties.VariableNames = h1;
Temp_outputAUC_Table = outputAUC_Table(1:end-1,:);
All_outputAUC_Table = outputAUC_Table(end,:);
Temp_outputAUC_Table = sortrows(Temp_outputAUC_Table,"Peptide");
outputAUC_Table = [Temp_outputAUC_Table;All_outputAUC_Table];
writetable(outputAUC_Table,outfilename,'Sheet','NormAUC');

%
outputstructure=cell(length(header),1); %Prep for Upset plot
for i = 1: length(pep)-1
    idx = logical(Count2(i,2:end));
    tempheader = header(2:end);
    tempheader(~idx) = ''; %deletes all structures not found on peptide
    emptcell=cell(length(header)-length(tempheader),1);
    tempheader=[tempheader,emptcell'];
    outputstructure=[outputstructure,tempheader'];
end
outputstructure=outputstructure(1:end-1,:);
outputstructure=[outputstructure,header(2:end)'];
outputstructure=cell2table(outputstructure(:,2:end));
t1=string(pep');
outputstructure.Properties.VariableNames = t1;
writetable(outputstructure,outfilename,'Sheet','UpSet'); %Upset table with 4 peptides as column names and list of present glycans as rows



% make plots
aa=figure(1);
aa.Position = [0 0 1920 1080];
glyClasses1=categorical(classes);
glyClasses1 = reordercats(glyClasses1,{'None','high_mannose','pauci','hybrid','complex','non-Fuc',...
    'coreFuc','termFuc','bothFuc','Neu5AcTerm','GalTerm','GlcNAcTerm','mono','bi','tri','tetra',...
    'other','bisect','non-bisect','lacdinac','lewis','slewis'});
bar(glyClasses1,Norm_AUC_classify') %Norm AUC, columns consist of Pep alone + 21 classes
ylim([0 100])
xlabel('Glycan Classes')
ylabel('Percent of AUC Class (%)')
title('Norm AUC Classify')
legend(pep,'Location','eastoutside')
saveas(aa,char(append(OutputFolder,'\',infilename,'_Ngly_Norm_AUC_Classify.png')))
set(gcf, 'Visible', 'off');

bb=figure(2);
bb.Position = [0 0 1920 1080];
classes(19)=[]; %remove non-bisecting
glyClasses2=categorical(classes);
glyClasses2 = reordercats(glyClasses2,{'None','high_mannose','pauci','hybrid','complex','non-Fuc',...
    'coreFuc','termFuc','bothFuc','Neu5AcTerm','GalTerm','GlcNAcTerm','mono','bi','tri','tetra',...
    'other','bisect','lacdinac','lewis','slewis'});
bar(glyClasses2,AUC_classify')
xlabel('Glycan Classes')
ylabel('Total AUC')
title('AUC Classify')
legend(pep)
saveas(bb,char(append(OutputFolder,'\',infilename,'_Ngly_AUC_Classify.png')))
set(gcf, 'Visible', 'off');

cc=figure(3);
cc.Position = [0 0 1920 1080];
bar(glyClasses2,Count_classify')
xlabel('Glycan Classes')
ylabel('# of HCD Scans')
title('Count Classify')
legend(pep)
saveas(cc,char(append(OutputFolder,'\',infilename,'_Ngly_Count_Classify.png')))
set(gcf, 'Visible', 'off');

%{
code for checking
# 1 To find peptide corresponding to HCD scan in Scan2
for i=1:length(Scan2{3,1})
pData.SGP{find(pData.Scan==Scan2{3,1}(i))}
end

# 2 To find glycopeptide corresponding to HCD scan in Scan2
for i=1:length(Scan2{3,4})
gpData.SGP{find(gpData.Scan_HCD==Scan2{3,4}(i))}
end 

# 3 To sum AUC for a set of glycopeptide HCD scans
AUCC
for i=1:length(Scan2{3,4})
AUCC=AUCC+gpData.AUC(find(gpData.Scan_HCD==Scan2{3,4}(i)))
end 

# 4 Find all scans corresponding to LeX
write scan number, AUC values for each
idx1 = (classifier(19,:) == 1)  % 19 us or Lex
idx2 = cellfun(@(x) ~isempty(x), (Scan2(1,:)))
idx3 = find(and(idx1, idx2))
AUC2{1,:}(~idx)=0

%}

a=1

