function [keepwhichprotein,allscoredist] = analyzepresearchdata(input,triggerdata,...
    allpepdata_pres,scoreoptions,rule)
% ANALYZEPRESEARCHDATA: analyze the result from the presearch and pick out
%     proteins that are more likely be fragmented in the experiment.
%
% Syntax:
% [keepwhichprotein,allscoredist] = analyzepresearchdata(input,triggerdata,...
%     allpepdata_pres,scoreoptions,rule)
%
% Input:
% input: structure.  Information necessary for scoring. See SCOREALLSPECTRA
%     for detail.
% triggerdata: structure. Scan number association in HCD trigger
%     experiment. See SEARCHSPECTRA for detail.
% allpepdata_pres: structure. The information in digested glycopeptide
%     file. See DIGESTFILEANALY for details. Note that field "varptm" has
%     been replaced by glycans used in pre-search.
% scoreoptions: structure. Settings and parameters for scoring. See
%     SEARCHSPECTRA and SCOREALLSPECTRA for detail.
% rule: double. The rule for identifying proteins to be searched in the
%     main search. Rules are:
%         1: proteins with peptides that have higher CombiES will be
%             selected first.
%         2: similar to rule 1, but the similarities between proteins are
%             considered. If the protein under consideration have no less
%             than half of its peptides identical to a protein that is
%             already in the list, this protein will be included, but not
%             counted.  Therefore int the final list the number of proteins
%             may exceed the limit.
%         3: for each protein the CombiES of all peptides that are above
%             0.6 are accumulated. The protein with the highest sum of
%             scores are chosen.
%
% Output:
% keepwhichprotein: n x 1 logical. The indices of proteins to be analyzed
%     in the main search. n equals to the number of proteins in the input.
% allscoredist: n x 10 numerical array. The distribution of glycopeptides'
%    CombiES of each protein. Each position represents a score range of
%    0.1. For example, allscoredist(100,7) means how many glycopeptides of
%    the 100th protein have CombiES greater than 0.6 and no higher than
%    0.7.
%
% Note:
% By default rule 3 is chosen in the program.
% For rule 1 and 2, for a protein to be considered, its must have at least
%     1 peptide whose CombiES >= 0.6, or at least 2 peptides whose
%     CombiES >= 0.6.
%
% Example:
% N/A
%
% Children function:
% N/A
%
% See Also:
% GROUPTRIGGEREDRESULT CALCULATECOMBIES

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 04/01/2020

topnprot_method5 = GlycoPATConstant.TopNProt_Method5_analyzepresearchdata;
CombiES_threshold = GlycoPATConstant.EScutoff_analyzepresearchdata;
Pepcov_threshold = GlycoPATConstant.Pepcovcutoff_analyzepresearchdata;

if ~input.dohcdtrigger
    errordlg('Presearch is available for HCD triggered experiment only.',...
        'Check input.');
    return
end

wtbar_analyzepresearchdata = waitbar(0,'Loading presearch results...');
keeptopnprot = input.presearchtopprotnum;  % Upper limit of the final list's size
numprotperbatch = input.doparacomp * input.numprotperworker;
% In pre-search, result files' name follow a certain pattern, which can be
%     calculated using algorithms. The detail of this algorithm can be
%     found in SAVESCORERESULTS
totalprotnum = length(allpepdata_pres.glypep);
keepwhichprotein = false(totalprotnum,1);
[~,f,e] = fileparts(input.outputfname);
totalfilenum = ceil(totalprotnum/numprotperbatch);
tgtfolderdir = dir(input.outputdir);
existingdatafilenames = {tgtfolderdir.name};
existingdatafilenames = existingdatafilenames(~[tgtfolderdir.isdir]);
presearchfilenames = {};
for ii = 1:totalfilenum
    protindse = [(ii-1)*numprotperbatch+1,...
        min(ii*numprotperbatch,totalprotnum)];
    tempoutputfilename = [f,'_presearch','_',num2str(protindse(1)),'_',num2str(protindse(2)),e];
    if ismember(tempoutputfilename,existingdatafilenames)
        presearchfilenames = [presearchfilenames;tempoutputfilename];
    else
        tempoutputfilename_regexppattern = [f,'_presearch','_',num2str(protindse(1)),'_',...
            num2str(protindse(2)),'_Part_[0-9]+of[0-9]+',e];
        fileexists = ~cellfun(@isempty,regexp(existingdatafilenames,tempoutputfilename_regexppattern));
        temppresearchfilenames = existingdatafilenames(fileexists);
        presearchfilenames = [presearchfilenames;temppresearchfilenames(:)];
    end
end

%% CALCULATE SCORE FOR EACH PEPTIDE SPECTRUM MATCH
switch rule
    case {1 2 3}
        allscoredist = zeros(totalprotnum,10);
        lower = 0:.1:.9;
        higher = .1:.1:1;
        for ii = 1:length(presearchfilenames)
            temppresearchfilename = fullfile(input.outputdir,presearchfilenames{ii});
            if exist(temppresearchfilename,'file')
                load(fullfile(input.outputdir,presearchfilenames{ii}),'result');
            else
                errordlg(['Pre-search file ',temppresearchfilename,' not found'],...
                    'Check data integrity');
                return
            end
            % Borrowed components from the result browser to calculate CombiES
            displaydataind = grouptriggeredresult(triggerdata.SASSO,...
                {result.ProteinID},[result.Scan],'regular');
            [CombiES,~] = calculateCombiES(result,displaydataind,...
                triggerdata.colnames);
            % Calculate the distribution of glycopeptide's CombiES
            for jj = 1:length(CombiES)
                representdispind = find(displaydataind(jj,:),1,'first');
                representdispind = displaydataind(jj,representdispind);
                thisprotid = str2num(result(representdispind).ProteinID);
                thisprotid = thisprotid(1);
                horipos = find(CombiES(jj) > lower & CombiES(jj) <= higher);
                if isempty(horipos)  % in this case CombiES = 0
                    horipos = 1;
                end
                allscoredist(thisprotid,horipos) = allscoredist(thisprotid,horipos) + 1;
            end
            waitbar(ii/totalfilenum * 0.9,wtbar_analyzepresearchdata);
        end
        waitbar(ii/totalfilenum * 0.9,wtbar_analyzepresearchdata,...
            'Analyzing presearch results...');
    case 4
        proteinscore = zeros(totalprotnum,1);
        for ii = 1:length(presearchfilenames)
            temppresearchfilename = fullfile(input.outputdir,presearchfilenames{ii});
            if exist(temppresearchfilename,'file')
                load(fullfile(input.outputdir,presearchfilenames{ii}),'result');
            else
                errordlg(['Pre-search file ',temppresearchfilename,' not found'],...
                    'Check data integrity');
                return
            end
            result = result([result.Fracpepfragmatch] >= 0.2);
            for jj = 1:length(result)
                tempprotid = str2num(result(jj).ProteinID);
                proteinscore(tempprotid(1)) = proteinscore(tempprotid(1)) + 1;
            end
            waitbar(ii/totalfilenum * 0.9,wtbar_analyzepresearchdata);
        end
        waitbar(ii/totalfilenum * 0.9,wtbar_analyzepresearchdata,...
            'Analyzing presearch results...');
    case 5
        proteinscore = zeros(totalprotnum,1);
        presearchresult = [];
        for ii = 1:length(presearchfilenames)
            temppresearchfilename = fullfile(input.outputdir,presearchfilenames{ii});
            if exist(temppresearchfilename,'file')
                load(fullfile(input.outputdir,presearchfilenames{ii}),'result');
                presearchresult = [presearchresult,result];
            else
                errordlg(['Pre-search file ',temppresearchfilename,' not found'],...
                    'Check data integrity');
                return
            end
            waitbar(ii/totalfilenum * 0.9,wtbar_analyzepresearchdata);
        end
        % Borrowed components from the result browser to calculate CombiES
        displaydataind = grouptriggeredresult(triggerdata.SASSO,...
            {presearchresult.ProteinID},[presearchresult.Scan],'regular');
        [CombiES,~] = calculateCombiES(presearchresult,displaydataind,...
            triggerdata.colnames);
        presearch_keepind = false(size(displaydataind,1),1);
        presearch_scans = [presearchresult(displaydataind(:,1)).Scan];
        presearch_prots = {presearchresult(displaydataind(:,1)).Protein};
        
        [~,~,ind2] = unique(presearch_scans);
        for ii = 1:max(ind2)
            temprowind = find(ind2 == ii);
            tempCombiES = CombiES(temprowind);
            [unitempCombiES,ind,~] = unique(tempCombiES);
            tempCombiES_sorted = flip(unitempCombiES);
            tempmaxPepcov = zeros(size(temprowind));
            for jj = 1:length(temprowind)
                tempmaxPepcov(jj) = max([presearchresult(displaydataind(temprowind(jj),:)).Pepcov]);
            end
            tempmaxPepcov_sorted = tempmaxPepcov(flip(ind));
            tempprots = presearch_prots(temprowind);
            tempprots_sorted = tempprots(flip(ind));
            protchosen = {};
            while ~isempty(tempCombiES_sorted) && (length(protchosen) < topnprot_method5)
                if (tempCombiES_sorted(1) >= CombiES_threshold) && (tempmaxPepcov_sorted(1) >= Pepcov_threshold)
                    presearch_keepind(temprowind(tempCombiES == tempCombiES_sorted(1))) = true;
                    protchosen = [protchosen;tempprots_sorted{1}];
                    protchosen = unique(protchosen);
                end
                tempCombiES_sorted(1) = [];
                tempmaxPepcov_sorted(1) = [];
                tempprots_sorted(1) = [];
            end
%             for jj = 1:min(topnprot_method5,length(tempCombiES_sorted))
%                 thistempCombiES_sorted = tempCombiES_sorted(jj);
%                 thistemp
%                 if thistempCombiES_sorted >= CombiES_threshold
%                     presearch_keepind(temprowind(tempCombiES == thistempCombiES_sorted)) = true;
%                 end
%             end
        end
        displaydataind_filtered = displaydataind(presearch_keepind,:);
        CombiES_filtered = CombiES(presearch_keepind);
        presearchdatasample_filtered = presearchresult(displaydataind_filtered(:,1));
        protid_filtered = zeros(length(presearchdatasample_filtered),1);
        for ii = 1:length(presearchdatasample_filtered)
            tempproteinid = str2num(presearchdatasample_filtered(ii).ProteinID);
            protid_filtered(ii) = tempproteinid(1);
        end
        [uniprotid,~,uniprotind2] = unique(protid_filtered);
        for ii = 1:max(uniprotind2)
            proteinscore(uniprotid(ii)) = sum(CombiES_filtered(uniprotind2 == ii))/...
                length(allpepdata_pres.glypep{uniprotid(ii)});
        end
    otherwise
        errordlg('Presearch rule not available.',...
            'Check input.');
end

%% RANK EACH PSM
switch rule
    case 1
        searchlowerrange = true;  % if CombiES >= 0.7 does not give enough proteins search
        % CombiES >= 0.6
        for ii = 10:-1:8
            currentgoodprotnum = sum(keepwhichprotein);  % Current number of selected protein
            tempallgoodprot = keepwhichprotein;
            goodprotind = allscoredist(:,ii) > 0;
            tempallgoodprot(goodprotind) = true;  % How many more proteins
            if sum(tempallgoodprot) > keeptopnprot  % If including proteins at this level will exceed limit
                [~,ind] = sort(allscoredist(:,ii),'descend');
                for jj = 1:length(ind)  % Add these proteins one by one. Proteins with more glycopeptides at
                    % this level first
                    if (currentgoodprotnum < keeptopnprot) && (~keepwhichprotein(ind(jj)))
                        keepwhichprotein(ind(jj)) = true;
                        currentgoodprotnum = currentgoodprotnum + 1;
                    end
                end
                searchlowerrange = false;
            else  % Including proteins at this level won't exceed limit
                keepwhichprotein(goodprotind) = true;  % Add them to list directly
            end
        end
        if searchlowerrange  % Search CombiES >= 0.6, the rest of the process is similar
            currentgoodprotnum = sum(keepwhichprotein);
            tempallgoodprot = keepwhichprotein;
            goodprotind = allscoredist(:,7) > 1;
            tempallgoodprot(goodprotind) = true;
            if sum(tempallgoodprot) > keeptopnprot
                [~,ind] = sort(allscoredist(:,7),'descend');
                for jj = 1:length(ind)
                    if (currentgoodprotnum < keeptopnprot) && (~keepwhichprotein(ind(jj)))
                        keepwhichprotein(ind(jj)) = true;
                        currentgoodprotnum = currentgoodprotnum + 1;
                    end
                end
            else
                keepwhichprotein(goodprotind) = true;
            end
        end
    case 2  % Method is similar to rule 1, but before adding a protein to the list,
        % similarity of this protein with those already in the list will be calculated.
        currentgoodprotnum = 0;
        for ii = 10:-1:8
            if currentgoodprotnum < keeptopnprot
                [~,ind] = sort(allscoredist(:,ii),'descend');
                for jj = 1:length(ind)
                    if (currentgoodprotnum < keeptopnprot) && (~keepwhichprotein(ind(jj))) &&...
                            allscoredist(ind(jj),ii) > 0
                        pickedglypep = allpepdata_pres.glypep(keepwhichprotein);
                        keepwhichprotein(ind(jj)) = true;
                        thisprotglypep = allpepdata_pres.glypep{ind(jj)};
                        partexist = 0;
                        % "Similarity" between protein A and B is defined by the proportion
                        %     of peptides of A whom can be found in B
                        for k = 1:length(pickedglypep)
                            partexist = max(partexist,sum(ismember(thisprotglypep,pickedglypep{k}))/...
                                length(thisprotglypep));
                        end
                        if partexist < 0.5  % Similarity >= 0.5 means the 2 proteins are similar
                            % therefore not counted
                            currentgoodprotnum = currentgoodprotnum + 1;
                        end
                    end
                end
            end
        end
        if currentgoodprotnum < keeptopnprot
            [~,ind] = sort(allscoredist(:,7),'descend');
            for jj = 1:length(ind)
                if (currentgoodprotnum < keeptopnprot) && (~keepwhichprotein(ind(jj))) &&...
                        allscoredist(ind(jj),7) > 1
                    pickedglypep = allpepdata_pres.glypep(keepwhichprotein);
                    keepwhichprotein(ind(jj)) = true;
                    thisprotglypep = allpepdata_pres.glypep{ind(jj)};
                    partexist = 0;
                    for k = 1:length(pickedglypep)
                        partexist = max(partexist,sum(ismember(thisprotglypep,pickedglypep{k}))/...
                            length(thisprotglypep));
                    end
                    if partexist < 0.5
                        currentgoodprotnum = currentgoodprotnum + 1;
                    end
                end
            end
        end
    case 3  % CombiES
        presearchscore = zeros(size(allscoredist,1),1);
        for ii = 1:size(allscoredist,1)
            if any(allscoredist(ii,9:10)) || allscoredist(ii,8) > 1
                presearchscore(ii) = allscoredist(ii,8:10) * [.7 .8 .9]';
            end
        end
        [~,ind] = sort(presearchscore,'descend');
        keepwhichprotein(ind(1:keeptopnprot)) = true;
    case {4,5}  % Pure PSM count | sum(CombiES)/length(glypeplib)
        [proteinscore_sorted,proteinscoreind] = sort(proteinscore,'descend');
        if any(proteinscore_sorted)
            keepwhichprotein(proteinscoreind(1:min(sum(proteinscore_sorted > 0),keeptopnprot))) = true;
        end
        allscoredist = proteinscore;
end
close(wtbar_analyzepresearchdata);
end