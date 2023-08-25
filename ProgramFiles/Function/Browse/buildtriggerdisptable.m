function [disptable,colnames,disptableind] = buildtriggerdisptable(scoredata,...
    SASSO,SASSO_method,disptableind)
% BUILDTRIGGERDISPTABLE: for trigger experiments, based on scores and scan
%     number association, find out where to put the results in display table.
%
% Syntax:
% [disptable,colnames,disptableind] = buildtriggerdisptable(scoredata,...
%     SASSO,SASSO_method,disptableind)
%
% Input:
% scoredata: result data, must be a table.
% SASSO: n x 5 numerical array, subsequent scan association data.
%     First column is trigger scan, after that are subsequent scans, the
%     last 2 columns are parent MS1 scan number and its retention time.
% SASSO_method: 1 x n cell array of strings, fragmentation methods of scans
%     in SASSO.
%     Check document of GROUPTRIGGEREDRESULT for another description.
% disptableind: n x m numerical array or empty, this variable tells where
%     to place each line of the result in the display table. If empty, this
%     function will build one. If not empty, the display table data will be
%     built using this variable.
%
% Output:
% disptable: m x n cell array, the display table data generated after
%     grouping with SASSO data.
% colnames: the name of each column in disptable. Including common ones and
%     frag. method specific ones. Common ones are: "Protein", "Expt", "Theo",
%     "Chg", "SGP", "combiES", "AUC", "Retime". Specific ones are: "Scan_X",
%     "peakLag_X", "P_X", "Top10_X", "%IonMatch_X" and "ES_X"
%     (X is fragmentation method name).
% disptableind: similar to the input with the same name, this variable
%     tells where each line of the result is placed in the output display table.
%
% Note:
% This function handles trigger experiment data only.
%
% Example:
% Set break points in BROWSEGUI
%
% Children function:
% GROUPTRIGGEREDRESULT
%

protein = scoredata.Protein;
scans = scoredata.Scan;
sgps = scoredata.SGP;
expt = scoredata.Mono;
mono = scoredata.Theo;
charge = scoredata.Charge;
peaklag = scoredata.PeakLag;
pvalue = scoredata.Pvalue;
top10 = scoredata.Top10;
percentionmatch = scoredata.PercentIonMatch;
enscore = scoredata.Enscore;
decoyes = scoredata.DecoyEnscore;
quant = scoredata.Quant;
nummethods = length(SASSO_method);
if isempty(disptableind)
    disptableind = grouptriggeredresult(SASSO,sgps,scans);
    colnames = {'Protein','Expt','Theo','Chg','SGP','CombiES','decoyCombiES','AUC','Retime'};
    for i = 1:nummethods
        tempmethod = SASSO_method{i};
        colnames = [colnames,['Scan','_',tempmethod],['peakLag','_',tempmethod],...
            ['P','_',tempmethod],['Top10','_',tempmethod],['%IonMatch','_',tempmethod],...
            ['ES','_',tempmethod],['decoyES','_',tempmethod]];
    end
else
    colnames = {};
end

disptbht = size(disptableind,1);
disptable = cell(disptbht,9);
for i = 1:disptbht  % "newtableind" contains zeroes in case scan got no hit
    tempind = disptableind(i,:);
    tempind = tempind(tempind > 0);
    disptable(i,1:5) = {protein{tempind(1)},expt(tempind(1)),...
        mono(tempind(1)),charge(tempind(1)),...
        sgps{tempind(1)}};
end
triggertable = cell(disptbht,7*nummethods);
triggeres = zeros(disptbht,nummethods);
triggerdecoyes = zeros(disptbht,nummethods);
triggerquant = zeros(disptbht,1);
triggerretime = zeros(disptbht,1);
for i = 1:disptbht
    tempnewtableind = disptableind(i,:);
    for j = 1:nummethods
        if tempnewtableind(j) > 0
            triggertable{i,(j-1)*7+1} = scans(tempnewtableind(j));
            triggertable{i,(j-1)*7+2} = peaklag(tempnewtableind(j));
            triggertable{i,(j-1)*7+3} = pvalue(tempnewtableind(j));
            triggertable{i,(j-1)*7+4} = top10(tempnewtableind(j));
            triggertable{i,(j-1)*7+5} = percentionmatch(tempnewtableind(j));
            triggertable{i,(j-1)*7+6} = enscore(tempnewtableind(j));
            triggertable{i,j*7} = decoyes(tempnewtableind(j));
            triggeres(i,j) = enscore(tempnewtableind(j));
            triggerdecoyes(i,j) = decoyes(tempnewtableind(j));
        end
    end
    whichfragmod = find(tempnewtableind,1);
    if ~isempty(whichfragmod)
        triggerquant(i) = quant(tempnewtableind(whichfragmod));
        tempretimes = SASSO(SASSO(:,whichfragmod) == ...
            triggertable{i,(whichfragmod-1)*7+1},end);
        triggerretime(i) = tempretimes(1);
    end
end
%% calculate combined ES
combiweight = zeros(size(SASSO_method));
for i = 1:length(SASSO_method)
    switch upper(SASSO_method{i})
        case 'CID'
            combiweight(i) = .25;
        case 'HCD'
            combiweight(i) = .35;
        case 'ETHCD'
            combiweight(i) = .4;
    end
end
combiES = zeros(size(disptableind,1),1);
decoycombiES = zeros(size(disptableind,1),1);
for i = 1:size(disptableind,1)
    tempdisptableind = disptableind(i,:);
    tempcombiweight = combiweight;
    tempcombiweight(tempdisptableind == 0) = 0;
    combiES(i) = triggeres(i,:)*tempcombiweight(:)/sum(tempcombiweight);
    decoycombiES(i) = triggerdecoyes(i,:)*tempcombiweight(:)/sum(tempcombiweight);
end



% combiES = triggeres*combiweight(:)/sum(combiweight);
% decoycombiES = triggerdecoyes*combiweight(:)/sum(combiweight);
for i = 1:disptbht
    disptable{i,6} = combiES(i);
    disptable{i,7} = decoycombiES(i);
    disptable{i,8} = triggerquant(i);
    disptable{i,9} = triggerretime(i);
end
disptable = [disptable,triggertable];
[~,ind] = sort(combiES,'descend');
disptable = disptable(ind,:);
disptableind = disptableind(ind,:);
end