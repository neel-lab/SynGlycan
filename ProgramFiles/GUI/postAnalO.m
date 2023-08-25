function postAnalO(pepList,glyList,pData,gpData,OutputFolder,Out_Loc)
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

% Get File Parts
[infilepath,infilename,infileext] = fileparts(Out_Loc);
OutputFolder=char(OutputFolder);
outfilename = strcat(OutputFolder,'\0_',infilename,'.xlsx');
outfilename = strrep(outfilename,'.xlsx','_OGlyPostAnal.xlsx');

% step 1: Find the average AUC for a given HCD scan in gp data and equally divide
% among hits from given HCD scan
gpHCD = gpData.Scan_HCD;
[~, ~, c] = unique(gpHCD);
gpQuant = gpData.AUC;
for i = 1:max(c)
    idx = find(c==i);
    gpData.AUC(idx)=mean(gpQuant(idx))/length(idx);
end

% targetGP='TGGTSTVDLGTSGTQAR';
% step 2. Find all glycopeptides of interest
[gpDataP,~,~]=cellfun(@(x) breakGlyPep(x), gpData.SGP, 'UniformOutput', false);
all_gp={};
for i=length(gpDataP):-1:1
    all_gp=[all_gp,gpDataP{i}.pep];
end
targetGP=unique(all_gp);

% step 3. Next for each target glycopeptide
pDataP = pData.SGP;
finalSiteCell={};
finalGlycanCell={};
for ii = 1:length(targetGP)
    temp_gpData = gpData;
    temp_pData = pData;

    for i=length(gpDataP):-1:1
        if ~strcmp(gpDataP{i,1}.pep,targetGP{ii}) %replaced contains with strcmp
            temp_gpData(i,:)=[];       % get rid of all except targetGP
        end
    end

    for i=length(pDataP):-1:1
        if ~strcmp(pDataP{i},targetGP{ii}) %replaced contains with strcmp
            temp_pData(i,:)=[];        % get rid of all except targetGP
        end
    end
    if ~isempty(temp_pData)
        pepCount = size(temp_pData,1);
        pepAUC = sum(temp_pData.Quant);
    else
        pepCount = 0;
        pepAUC = 0;
    end

    AUC=zeros(length(targetGP{ii}),size(glyList,1));     % to store AUC for each peptide/glycopeptide
    Count=zeros(length(targetGP{ii}),size(glyList,1));   % to count number of spectra for each (glyco)peptide

    % step 4: Summarize the data for O-glycans based on structure and position
    % of occupancy
    for i = 1:size(temp_gpData)
        [~,gly,~]= breakGlyPep(temp_gpData.SGP{i});
        for j=1:size(gly,2)
            glycan_idx=find(strcmp(glyList,gly(j).struct));
            pos_idx=gly(j).pos;
            AUC(pos_idx,glycan_idx)=AUC(pos_idx,glycan_idx)+temp_gpData.AUC(i);
            Count(pos_idx,glycan_idx)=Count(pos_idx,glycan_idx)+1;
        end
    end

    % step 5: make plots
    SiteSum=sum(AUC,2);      % sum rows
    GlycanSum=sum(AUC,1);    % sum Columns
    GlycanSum=[pepAUC, GlycanSum];
    glyList_all=[' ';glyList];
    dd = figure(4);
    dd.Position = [0 0 1920 1080];
    bar(SiteSum)
    X_label=cell(length(targetGP{ii}),1);
    set(gca, 'XTick', 1:length(targetGP{ii}))
    for i =1:length(targetGP{ii})
        X_label{i}=targetGP{ii}(i);
    end
    set(gca,'xticklabel',X_label)
    xlabel('Site Number on targetGP')
    ylabel('AUC')
    title('AUC Per Peptide Site')
    saveas(dd,strcat(char(OutputFolder),'\',infilename,'_',targetGP{ii},'_OGP_AUC_SiteSum.png'))
    set(gcf, 'Visible', 'off');

    ee = figure(5);
    ee.Position = [0 0 1920 1080];
    bar(GlycanSum)
    set(gca, 'XTick', 1:length(glyList_all))
    set(gca,'xticklabel', glyList_all)
    title('AUC Per Glycan Structure')
    xlabel('O-Glycan Structure')
    ylabel('AUC')
    title('AUC Per O-Glycan Structure')
    saveas(ee,strcat(char(OutputFolder),'\',infilename,'_',targetGP{ii},'_OGP_AUC_GlycanSum.png'))
    set(gcf, 'Visible', 'off');

    SiteCell=cell(2,length(targetGP{ii}));
    for i=1:length(targetGP{ii})
        SiteCell(1,i)=cellstr(targetGP{ii}(i));
        SiteCell(2,i)=num2cell(SiteSum(i));
    end
    GlycanCell=cell(2,length(targetGP));
    for i=1:size(glyList_all,1)
        GlycanCell(1,i)=glyList_all(i);
        GlycanCell(2,i)=num2cell(GlycanSum(i));
    end

    size1=size(finalSiteCell);
    size2=size(SiteCell);
    % Pad the smaller cell array with empty cells to make their sizes compatible
    if size1(2) < size2(2)
        finalSiteCell = [finalSiteCell, repmat({''}, size1(1), size2(2) - size1(2))];
    elseif size1(2) > size2(2)
        SiteCell = [SiteCell, repmat({''}, size2(1), size1(2) - size2(2))];
    end
    % Concatenate the two cell arrays vertically
    finalSiteCell = [finalSiteCell; SiteCell];


    size1=size(finalGlycanCell);
    size2=size(GlycanCell);
    % Pad the smaller cell array with empty cells to make their sizes compatible
    if size1(2) < size2(2)
        finalGlycanCell = [finalGlycanCell, repmat({''}, size1(1), size2(2) - size1(2))];
    elseif size1(2) > size2(2)
        GlycanCell = [GlycanCell, repmat({''}, size2(1), size1(2) - size2(2))];
    end
    % Concatenate the two cell arrays vertically
    finalGlycanCell = [finalGlycanCell; GlycanCell];
end


% step 4: Consider data for peptide with no glycan

writecell(finalSiteCell,char(outfilename),'Sheet','AUC_for_each_position');
writecell(finalGlycanCell,char(outfilename),'Sheet','AUC_for_each_glycan');

end
