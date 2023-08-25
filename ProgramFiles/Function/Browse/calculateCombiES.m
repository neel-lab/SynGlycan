function [CombiES,CombiDecoyES] = calculateCombiES(scoredata,displaydataind,...
    fragmethods)
% CALCULATECOMBIES: calculate combiES and combiDecoyES for result of HCD
%     triggered experiment.
% 
% Syntax:
% [CombiES,CombiDecoyES] = calculateCombiES(scoredata,displaydataind,fragmethods,workingmode)
% 
% Input:
% 
% 
% Output:
% 
% 
% Note:
% 
%
% Example:
% 
%
% Children function: 
% 

% GlycoPAT 2 authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2020, Research Foundation for State University of New York. All rights reserved

alles = [scoredata.Enscore];
alldecoyes = [scoredata.DecoyEnscore];
numfragmethods = length(fragmethods);
combiweight = zeros(numfragmethods,1);
for ii = 1:numfragmethods
    switch upper(fragmethods{ii})
        case 'CID'
            combiweight(ii) = GlycoPATConstant.combiweight.CID;
        case 'HCD'
            combiweight(ii) = GlycoPATConstant.combiweight.HCD;
        case 'ETHCD'
            combiweight(ii) = GlycoPATConstant.combiweight.ETHCD;
    end
end
if iscell(displaydataind)  % meaning isomers are grouped
    CombiES = cell(size(displaydataind,1),1);
    CombiDecoyES = cell(size(displaydataind,1),1);
    for ii = 1:size(displaydataind,1)
        numisomers = length(displaydataind{ii,1});
        tempes = zeros(numisomers,numfragmethods);
        tempdecoyes = zeros(numisomers,numfragmethods);
        tempglypepind = zeros(numisomers,numfragmethods);
        tempCombiES = zeros(numisomers,1);
        tempCombiDecoyES = zeros(numisomers,1);
        for jj = 1:numisomers
            for kk = 1:numfragmethods
                tempglypepind2 = displaydataind{ii,kk}(jj);
                if tempglypepind2 > 0
                    tempes(jj,kk) = alles(tempglypepind2);
                    tempdecoyes(jj,kk) = alldecoyes(tempglypepind2);
                    tempglypepind(jj,kk) = tempglypepind2;
                end
            end
            tempCombiES(jj) = sum(tempes(jj,:) * combiweight(:)) / ...
                sum(combiweight(tempglypepind(jj,:) > 0));
            tempCombiDecoyES(jj) = sum(tempdecoyes(jj,:) * combiweight(:)) / ...
                sum(combiweight(tempglypepind(jj,:) > 0));
        end
        CombiES{ii} = tempCombiES;
        CombiDecoyES{ii} = tempCombiDecoyES;   
    end    
elseif isnumeric(displaydataind)  % meaning each isomer occupy a line
    CombiES = zeros(size(displaydataind,1),1);
    CombiDecoyES = zeros(size(displaydataind,1),1);
    for ii = 1:size(displaydataind,1)
        tempglypepind = displaydataind(ii,:);
        tempes = zeros(1,numfragmethods);
        tempdecoyes = zeros(1,numfragmethods);
        for jj = 1:numfragmethods
            if tempglypepind(jj) > 0
                tempes(jj) = alles(tempglypepind(jj));
                tempdecoyes(jj) = alldecoyes(tempglypepind(jj));
            end
        end
        CombiES(ii) = sum(tempes * combiweight(:)) / ...
            sum(combiweight(tempglypepind > 0));
        CombiDecoyES(ii) = sum(tempdecoyes * combiweight(:)) / ...
            sum(combiweight(tempglypepind > 0));
    end
end
end