function escutoff = fdrfilter(scoredata,displayoptions,displaydataind,...
    fdrfiltervalue,fdrfragmode)
% FDRFILTER: calculate ensemble score cutoffs at a given FDR level.
%
% Syntax:
% escutoff = fdrfilter(scoredata,displayoptions,displaydataind,...
%     fdrfiltervalue,fdrfragmode)
%
% Input:
% scoredata: 1 x n structure. Scoring data.
% displayoptions: structure. BrowseGUI display options. Currently include
%     fields: "scorefname" (result file name), "exptdatafname" (experiment
%     data file name), "combineisomer" (combine isomers as groups Y/N),
%     "scoreintdata" (a structure holding scoring parameters).
% displaydataind: m x n cell array of p x 1 double or m x n double. The
%     serial numbers of the results displayed in BrowseGUI table.
% fdrfiltervalue: double. FDR threshold.
% fdrfragmode: string. Working mdoe, currently accepted ones are:
%     "CombiES" (HCD triggered experiment only), "Every hit" (put all
%     identifications together to calculate), "Each frag. individually"
%     (each fragmentation method has a separate cutoff) and any
%     fragmentation method that appeared.
%
% Output:
% escutoff: 1 x n double. Calculated ensemble score cutoff for given FDR
%     level. n = 1 for most cases except "fdrfragmode" is "Each frag.
%     individually", where the corresponding fragmentation modes can be
%     found in "displayoptions.scoreintdata.analyzefragmode".
%
% Note:
% N/A
%
% Example:
% N/A
%
% Children function:
% N/A

% GlycoPAT authors: Kai Cheng, Sriram Neelamegham
%(c) 2020, Research Foundation for State University of New York. All rights reserved

scoreintdata = displayoptions.scoreintdata;
allfragmode = upper({scoredata.Fragmode});
alles = [scoredata.Enscore];
alldecoyes = [scoredata.DecoyEnscore];

if scoreintdata.scoreoptions.isHCDtrigger  % trigger
    if strcmpi(fdrfragmode,'COMBIES')
        SASSO_method = scoreintdata.colnames;
        numfragmethods = length(SASSO_method);
        combiweight = zeros(numfragmethods,1);
        for ii = 1:numfragmethods
            switch upper(SASSO_method{ii})
                case 'CID'
                    combiweight(ii) = GlycoPATConstant.combiweight.CID;
                case 'HCD'
                    combiweight(ii) = GlycoPATConstant.combiweight.HCD;
                case 'ETHCD'
                    combiweight(ii) = GlycoPATConstant.combiweight.ETHCD;
            end
        end
        if displayoptions.combineisomer
            allcombies = zeros(sum(cellfun(@length,displaydataind(:,1))),1);
            allcombidecoyes = zeros(size(allcombies));
            writeind = 1;
            for ii = 1:size(displaydataind,1)
                % triggered + combineisomer means "displaydataind" is n x m
                % cell array of p x 1 double. Each row is an isomer group,
                % each column is a frag method. Elements in "p" are
                % different isomers. For each isomer group, extract their
                % ES and dES, then calculate combiES and combidecoyES.
                tempes = zeros(length(displaydataind{ii,1}),numfragmethods);
                tempdecoyes = zeros(size(tempes));  % rows are isomers, columns are frag methods
                dataispresent = false(size(tempes));
                for jj = 1:numfragmethods
                    for kk = 1:length(displaydataind{ii,jj})
                        if displaydataind{ii,jj}(kk) > 0
                            tempes(kk,jj) = alles(displaydataind{ii,jj}(kk));
                            tempdecoyes(kk,jj) = alldecoyes(displaydataind{ii,jj}(kk));
                            dataispresent(kk,jj) = true;
                        end
                    end
                end
                tempcombies = zeros(length(displaydataind{ii,jj}),1);
                tempdecoycombies = zeros(length(displaydataind{ii,jj}),1);
                for jj = 1:length(displaydataind{ii,jj})
                    % If an identification is missing in a trigger group,
                    % combiES will use existing scores. But the weight of
                    % the score must be adjusted accordingly.
                    tempcombies(jj) = tempes(jj,:) * combiweight / ...
                        sum(combiweight(dataispresent(jj,:)));
                    tempdecoycombies(jj) = tempdecoyes(jj,:) * combiweight / ...
                        sum(combiweight(dataispresent(jj,:)));
                end
                allcombies(writeind:writeind + jj - 1) = tempcombies;
                allcombidecoyes(writeind:writeind + jj - 1) = tempdecoycombies;
                writeind = writeind + jj;
            end
        else  % do not combine isomers - "displaydataind" is matrix
            allcombies = zeros(size(displaydataind,1),1);
            allcombidecoyes = zeros(size(allcombies));
            for ii = 1:size(displaydataind,1)
                tempes = zeros(1,numfragmethods);
                tempdecoyes = zeros(size(tempes));  % rows are isomers, columns are frag methods
                for jj = 1:numfragmethods
                    if displaydataind(ii,jj) > 0
                        tempes(jj) = alles(displaydataind(ii,jj));
                        tempdecoyes(jj) = alldecoyes(displaydataind(ii,jj));
                    end
                end
                allcombies(ii) = tempes * combiweight / ...
                    sum(combiweight(displaydataind(ii,:) > 0));
                allcombidecoyes(ii) = tempdecoyes * combiweight / ...
                    sum(combiweight(displaydataind(ii,:) > 0));
            end
        end
    end
end

enscorecounter = zeros(101,1);
fdrenscorecounter = zeros(101,1);
if strcmpi(fdrfragmode,'Every hit')
    for ii = 1:101
        enscorecounter(ii) = sum(alles >= (ii-1)*.01);
        fdrenscorecounter(ii) = sum(alldecoyes >= (ii-1)*.01);
    end
elseif strcmpi(fdrfragmode,'Each frag. individually')
    fragmethods = scoreintdata.scoreoptions.analyzefragmode;
    enscorecounter = zeros(101,length(fragmethods));
    fdrenscorecounter = zeros(101,length(fragmethods));
    for ii = 1:length(fragmethods)
        goodfragmode = strcmpi(allfragmode,fragmethods{ii});
        tempes = alles(goodfragmode);
        tempdes = alldecoyes(goodfragmode);
        for jj = 1:101
            enscorecounter(jj,ii) = sum(tempes >= (jj-1)*.01);
            fdrenscorecounter(jj,ii) = sum(tempdes >= (jj-1)*.01);
        end
    end
elseif strcmpi(fdrfragmode,'CombiES')
    for ii = 1:101
        enscorecounter(ii) = sum(allcombies >= (ii-1)*.01);
        fdrenscorecounter(ii) = sum(allcombidecoyes >= (ii-1)*.01);
    end
else  % specified a fragmentation mode
    goodfragmode = strcmpi(allfragmode,fdrfragmode);
    tempes = alles(goodfragmode);
    tempdes = alldecoyes(goodfragmode);
    for ii = 1:101
        enscorecounter(ii) = sum(tempes >= (ii-1)*.01);
        fdrenscorecounter(ii) = sum(tempdes >= (ii-1)*.01);
    end
end
FDRratio = fdrenscorecounter ./ enscorecounter;
FDRratio(arrayfun(@isnan,FDRratio)) = 0;
badFDR = FDRratio > fdrfiltervalue;
% result from FIND still fails FDR cut off but + 1 will get you a good one
escutoff = ones(1,size(FDRratio,2));
for ii = 1:size(FDRratio,2)
    escutoff(ii) = (find(badFDR(:,ii),1,'last') + 1) / 100;
end
escutoff(escutoff > 1) = 1;
end