function [pmposstorage,pmindstorage,directionseq] = smartPMdispenser(pmsgroup,...
    bondmap,distance,directionseq,workingmode)
% SMARTPMDISPENSER: smart penpendicular monosac. structure dispenser. Helps
%     arrange perpendicular glycan sub-structures within the bone structure
%
% Syntax:
% [pmsposstorage,pmsindstorage,directionseq] = smartPMdispenser(pmsgroup,...
%     bondmap,distance,directionseq,workingmode)
%
% Input:
% pmsgroup: n x 2 cell array. 1st column contains the serial number of the
%     monosac. that have PM attached to it, 2nd column contains the serial
%     number of the corresponding PM's. Elements in 2nd column are m x 2
%     double. 1st column is the serial number of PM and all its children
%     monosac., 2nd is instruction for perpendicular orientation: 0 for auto,
%     1 for left (counterclockwise), 2 for right (clockwise).
% bondmap: m x m double. Linkage map of the glycan.
% distance: 1 x m double. Distance of the monosac. in the glycan.
% directionseq: 1 x m double. Orientation of the monosac. in the glycan.
% workingmode: double. 1 for automatic distribution, 2 for linkage based
%     distribution.
%
% Output:
% pmposstorage: n x 2 cell array of m x 2 double. Positioning info. of the
%     PMs. The 1st column contains the position info for the PMs to be
%     drawn upwards, 2nd column contains those to be drawn downwards.
% pmindstorage: n x 2 cell array of m x 1 double. Serial numbers of the PMs
%     to be drawn upwards or downwards. 1st column contains upward ones,
%     2nd column contains downwards ones.
% directionseq: 1 x m double. Orientation of the monosac. after
%     upwards/downwards adjustment.
%
% Note:
% In mode 1, if more than one PM is directly attached to the same monosac.
%     in the bone structure, SMARTPMDISPENSER will automatically dispense
%     the PM subtrees evenly around the root monosac.: the first points
%     down, second points up, third points down again.
%
% Example:
% N/A. Set breakpoints in main program to see how this function works.
%
% Children function:
% DRAWRAWTREE, BRANCHEQUALIZER
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2020, Research Foundation for State University of New York. All rights reserved
%


pmposstorage = cell(size(pmsgroup));
pmindstorage = cell(size(pmsgroup));
for i = 1:size(pmsgroup,1)  % one subtree at a time
    tmpsubtreeroot = pmsgroup{i,1};  % temp sub tree root
    subtreemsind = pmsgroup{i,2};
    allarms = [];
    uparms = [];
    downarms = [];
    defaultarms = [];
    for j = 1:length(subtreemsind)
        allarms = [allarms;subtreemsind{j}(:,1)];
        if subtreemsind{j}(1,2) == 1  % -PL
            uparms = [uparms;j];
        elseif subtreemsind{j}(1,2) == 2  % -PR
            downarms = [downarms;j];
        elseif subtreemsind{j}(1,2) == 0  % -P or default, distribute evenly
            defaultarms = [defaultarms;j];
        end
    end
    allbranches = find(bondmap(tmpsubtreeroot,:));  % all MS that connects to root
    mainbranchmsind = setdiff(allbranches,allarms);  % decide main branch,
    switch workingmode
        case 1
            % remaining default PMS will be distributed so that eventually num of
            % uparms will be equal or bigger than downarms only by 1
            if ~isempty(defaultarms)
                numuparms = length(uparms);
                numdownarms = length(downarms);
                numallarms = numuparms + numdownarms + length(defaultarms);
                uadiff = numallarms - numuparms*2;  % uadiff > 0 means need to put PMS to the upper arm
                dadiff = numallarms - numdownarms*2;
                if dadiff > 0
                    cut = min(ceil(dadiff/2),length(defaultarms));
                    downarms = [downarms;defaultarms(1:cut)];
                    defaultarms(1:cut) = [];
                end
                if uadiff > 0
                    uparms = [uparms;defaultarms];
                end
            end            
        case 2
            while ~isempty(defaultarms)
                currentarmmsind = subtreemsind{defaultarms(1)}(1);  % current branch
                if currentarmmsind >= mainbranchmsind  % clockwise
                    downarms = [defaultarms(1);downarms];
                else  % counterclockwise
                    uparms = [uparms;defaultarms(1)];
                end
                defaultarms(1) = [];
            end
    end
    if ~isempty(uparms)
        alluparms = subtreemsind(uparms);
        tmparmcont = [];
        for j = 1:numel(alluparms)
            tmparmcont = [tmparmcont;alluparms{j}(:,1)];
        end
        tmparmbondmap = bondmap([tmpsubtreeroot;tmparmcont],[tmpsubtreeroot;tmparmcont]);
        tmparmdistance = distance([tmpsubtreeroot;tmparmcont]);
        uparmmspos = drawrawtree(tmparmbondmap,tmparmdistance,[]);
        uparmmspos(:,1) = uparmmspos(:,1) - uparmmspos(1,1);
        uparmmspos(:,2) = uparmmspos(:,2) - uparmmspos(1,2);
        uparmmspos = branchequalizer(uparmmspos,tmparmbondmap);
        tmprow = -uparmmspos(:,1);
        uparmmspos(:,1) = uparmmspos(:,2);
        uparmmspos(:,2) = tmprow;
        uparmmspos = uparmmspos(2:end,:);
        pmposstorage{i,1} = uparmmspos;
        temppmsind = pmsgroup{i,2}(uparms);
        pmsind = [];
        for j = 1:numel(temppmsind)
            pmsind = [pmsind;temppmsind{j}(:,1)];
        end
        pmindstorage{i,1} = pmsind;
    end
    if ~isempty(downarms)
        alldownarms = subtreemsind(downarms);
        tmparmcont = [];
        for j = 1:numel(alldownarms)
            tmparmcont = [tmparmcont;alldownarms{j}(:,1)];
        end
        tmparmbondmap = bondmap([tmpsubtreeroot;tmparmcont],[tmpsubtreeroot;tmparmcont]);
        tmparmdistance = distance([tmpsubtreeroot;tmparmcont]);
        downarmmspos = drawrawtree(tmparmbondmap,tmparmdistance,[]);
        downarmmspos(:,1) = downarmmspos(:,1) - downarmmspos(1,1);
        downarmmspos(:,2) = downarmmspos(:,2) - downarmmspos(1,2);
        downarmmspos = branchequalizer(downarmmspos,tmparmbondmap);
        tmprow = downarmmspos(:,1);
        downarmmspos(:,1) = downarmmspos(:,2);
        downarmmspos(:,2) = tmprow;
        downarmmspos = downarmmspos(2:end,:);
        pmposstorage{i,2} = downarmmspos;  % turn PMS subtree 90 deg
        temppmsind = pmsgroup{i,2}(downarms);
        pmsind = [];
        for j = 1:numel(temppmsind)
            pmsind = [pmsind;temppmsind{j}(:,1)];
        end
        pmindstorage{i,2} = pmsind;
    end
end

end