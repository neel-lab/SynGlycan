function [reformed,msposfix] = branchreform(thisgly)
% BRANCHREFORM: Rearrange branches in glycan structure according to linkage
%     information.
%
% Syntax:
% [reformed,msposfix] = branchreform(thisgly)
%
% Input:
% thisgly: string. SGP2.0 sequence of glycan.
%
% Output:
% reformed: string. Condensed IUPAC sequence of glycan after rearrangement.
% msposfix: 1 x n double. The position shift of each monosaccharide. For
%     example, if a glycan "{n{s}{h{s}}" has been rearranged into
%     "{n{h{s}}{s}}, the first "n" stays, so msposfix(1) = 0, first "{s}" was
%     moved to the end, position moved from 2 to 4, so msposfix(2) = 2.
%     "{h{s}}" changed from 3rd and 4th to 2nd and 3rd, so msposfix(3) and
%     msposfix(4) both = -1.
%     If no rearrangement happens, this variable is [0,0,0,...].
%
% Note:
% Glycan structure is sorted by ascending order, first anomeric carbon
%     (alpha, beta,...) then hydroxyl group(1, 2, 3,...).
%
% Example:
% [reformed,msposfix]=branchreform('{GalNAc(??-?){Neu5Ac(a2-6)}{Gal(b1-3){Neu5Ac(a2-6)}}}')
% 
% reformed =
% 
%     '{GalNAc(??-?){Gal(b1-3){Neu5Ac(a2-6)}}{Neu5Ac(a2-6)}}'
% 
% 
% msposfix =
% 
%      0     2    -1    -1
%
% Children function:
% GLYTREETRACKER
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2020, Research Foundation for State University of New York. All rights reserved
%

%% Preprocessing - block user customized information
thisgly = strrep(thisgly,'[','');
thisgly = strrep(thisgly,']','');
thisglynocurl = thisgly;
[optvalstart,optvalend] = regexp(thisgly,'["''].*?["'']','start','end');
for ii = 1:length(optvalstart)
    thisglynocurl(optvalstart(ii):optvalend(ii)) = ['"',repmat('?',1,optvalend(ii)-optvalstart(ii)-1),'"'];
end
%% Standard monosaccharide distance calculation process
indtemp = strfind(thisglynocurl,'{');
levelindex = zeros(2,length(thisgly));
indtemp2 = strfind(thisglynocurl,'}');
if ~ismember(1,indtemp)
    error('Does glycan starts with a "{"?');
end
levelindex(1,indtemp) = 1;
levelindex(1,indtemp2) = -1;
levelindex(2,1) = 1;
for ii = 2:size(levelindex,2)
    levelindex(2,ii) = levelindex(2,ii-1) + levelindex(1,ii);
end
%% Read linkage info
wholestr = zeros(1,length(thisgly));
wholestr(indtemp) = 1;
wholestr(indtemp2) = 1;
allbond = cell(sum(wholestr)/2,1);
writeind = 1;
tempstr = '';
readind = 1;
while readind <= length(wholestr)
    if wholestr(readind) ~= 1  % gather character to form the monosac. string
        tempstr = [tempstr,thisgly(readind)];
    else
        if ~isempty(tempstr)
            bondinfo = regexp(tempstr,'[a-z?][\d?]+-[\d?]+','match');
            allbond{writeind} = bondinfo{1};
            tempstr = '';
            writeind = writeind + 1;
        end
    end
    readind = readind + 1;
end

%% Calculate distance of each monosac, build bondmap
letterindex = zeros(1,length(thisgly));
letterindex(regexp(thisglynocurl,'[^{}]')) = 1;
distance = letterindex.*levelindex(2,:);
distance = distance(indtemp+1);  % all monosac's
if length(indtemp) > 1
    for ii = 2:length(indtemp)
        letterindex(indtemp(ii):end) = letterindex(indtemp(ii):end)/(ii-1)*ii;
    end
end
bondmap = zeros(length(distance));
readind = 1;
while readind < length(distance)
    if (distance(readind + 1) - distance(readind)) == 1
        bondmap(readind,readind + 1) = 1;  % consecutive numbers indicate bond
        readind = readind + 1;
    elseif (distance(readind + 1) - distance(readind)) < 1  % if chain is broken, go back to find its fork point
        thisind = distance(readind + 1);  % where it's broken
        itsforkpt = find(distance(1:readind) == thisind - 1,1,'last');  % where is the fork point
        bondmap(itsforkpt,readind + 1) = 1;  % mark this bond
        readind = readind + 1;  % keep going on
    end
end
%% Identify fork position - this is where rearrangement happens
isfork = find(sum(bondmap,2) > 1);
msseq = 1:size(bondmap,1);
% For all downstream monosac. at the fork point, read their linkage
%     information, then sort.
if any(isfork)
    for ii = length(isfork):-1:1
        forkchildren = find(bondmap(isfork(ii),:));
        branches = cell(size(forkchildren));
        bonds = cell(length(forkchildren),1);
        branchpos = zeros(length(forkchildren),2);
        for jj = 1:length(forkchildren)
            branches{jj} = glytreetracker(bondmap,forkchildren(jj),[],'down');
            branchpos(jj,1) = min(find(letterindex == forkchildren(jj))) - 1;
            branchpos(jj,2) = find(levelindex(2,branchpos(jj,1):end) == levelindex(2,branchpos(jj,1))-1,1) + branchpos(jj,1)-1;
            tempbond = allbond{forkchildren(jj)};
            carbon = regexp(tempbond,'[a-z?]','match');
            nrnum = regexp(tempbond,'[0-9?]+','match');
            bonds{jj,1} = carbon{1};
            bonds{jj,2} = nrnum{end};
        end
        tempind = 1:size(bonds,1);
        branchseq = zeros(size(tempind));
        ind1 = ones(size(bonds,1),1);
        writeind = 1;
        for jj = 1:max(ind1)
            [~,~,ind2] = unique(bonds(ind1 == jj,2));
            indperletter = tempind(ind1 == jj);
            for l = 1:max(ind2)
                indpernumber = indperletter(ind2 == l);
                branchseq(writeind:writeind + length(indpernumber) - 1) = indpernumber;
                writeind = writeind + length(indpernumber);
            end
        end
        % Branch reshuffle - swapping glycan subtrees 
        % Glycan sequence has 3 parts: head, middle and tail. Swapping
        % happens in "middle" section. Combine head - middle (new) - tail
        % will give you the revised structure.
        newbranches = branches(branchseq);
        thisglyhead = thisgly(1:min(min(branchpos))-1);
        thisglytail = thisgly(max(max(branchpos))+1:end);
        thisglymiddle = thisgly(min(min(branchpos)):max(max(branchpos)));
        letterindhead = letterindex(1:min(min(branchpos))-1);
        letterindtail = letterindex(max(max(branchpos))+1:end);
        letterindmiddle = letterindex(min(min(branchpos)):max(max(branchpos)));
        lvlindhead = levelindex(:,1:min(min(branchpos))-1);
        lvlindtail = levelindex(:,max(max(branchpos))+1:end);
        lvlindmiddle = levelindex(:,min(min(branchpos)):max(max(branchpos)));
        msseqhead = msseq(1:min(cellfun(@min,newbranches))-1);
        msseqtail = msseq(max(cellfun(@max,newbranches))+1:end);
        msseqmiddle = [];
        for jj = 1:length(newbranches)
            msseqmiddle = [msseqmiddle,msseq(newbranches{jj})];
        end
        msseq = [msseqhead,msseqmiddle,msseqtail];
        minmsind = cellfun(@min,branches);
        maxmsind = cellfun(@max,branches);
        bondmaphead = bondmap(1:min(minmsind)-1,:);
        bondmaptail = bondmap(max(maxmsind)+1:end,:);
        branchpos = branchpos-min(min(branchpos))+1;
        newglymiddle = [];
        newletterindmiddle = [];
        newlvlindmiddle = [];
        newbondmapmiddle = [];
        newbranchlength = cellfun(@length,branches(branchseq));
        newforkchildren = forkchildren;
        for jj = 2:length(forkchildren)
            newforkchildren(jj) = newforkchildren(jj-1) + newbranchlength(jj-1);
        end
        bondmapforkline = zeros(1,size(bondmap,2));
        bondmapforkline(newforkchildren) = 1;
        bondmaphead(end,:) = bondmapforkline;
        bondmapfix = newforkchildren - forkchildren(branchseq);
        % The swapping takes place here
        for jj = 1:length(branchseq)
            newglymiddle = [newglymiddle,thisglymiddle(branchpos(branchseq(jj),1):branchpos(branchseq(jj),2))];
            newletterindmiddle = [newletterindmiddle,letterindmiddle(branchpos(branchseq(jj),1):branchpos(branchseq(jj),2))];
            newlvlindmiddle = [newlvlindmiddle,lvlindmiddle(:,branchpos(branchseq(jj),1):branchpos(branchseq(jj),2))];
            tempbondmapmiddle = bondmap(branches{branchseq(jj)},:);
            if bondmapfix(jj) < 0
                tempbondmapmiddle = [tempbondmapmiddle(:,-bondmapfix(jj)+1:end),zeros(size(tempbondmapmiddle,1),-bondmapfix(jj))];
            elseif bondmapfix(jj) > 0
                tempbondmapmiddle = [zeros(size(tempbondmapmiddle,1),bondmapfix(jj)),tempbondmapmiddle(:,1:size(tempbondmapmiddle,2)-bondmapfix(jj))];
            end
            newbondmapmiddle = [newbondmapmiddle;tempbondmapmiddle];
        end
        thisgly = [thisglyhead,newglymiddle,thisglytail];
        letterindex = [letterindhead,newletterindmiddle,letterindtail];
        levelindex = [lvlindhead,newlvlindmiddle,lvlindtail];
        bondmap = [bondmaphead;newbondmapmiddle;bondmaptail];
    end
end
reformed = thisgly;
[~,ind] = sort(msseq);
msposfix = ind - (1:length(msseq));
end