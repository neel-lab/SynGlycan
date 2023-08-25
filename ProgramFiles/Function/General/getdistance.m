function distance = getdistance(sgpseq)
% GETDISTANCE: for each monosaccharide, get the distance from glycan
% reducing end.
%
% Syntax:
% distance = getdistance(sgpseq)
%
% Input:
% sgpseq: glycan SGP 1.0 sequence.
%
% Output:
% distance: 1 x n numerical array, distance of each monosaccharide. n
% equals to number of monosaccharides in glycan.
%
% Note:
% N/A
%
% Example:
% distance = getdistance('{n{f}{h{s}}}')
%
% distance =
%
%      1     2     2     3
%
% Children function:
% N/A

thisgly = sgpseq;
indtemp = strfind(thisgly,'{');
levelindex = zeros(2,length(thisgly));
indtemp2 = strfind(thisgly,'}');
levelindex(1,indtemp) = 1;
levelindex(1,indtemp2) = -1;
levelindex(2,1) = 1;
for j = 2:size(levelindex,2)
    levelindex(2,j) = levelindex(2,j-1) + levelindex(1,j);
end
letterindex = zeros(1,length(thisgly));
letterindex(regexp(thisgly,'[^{}]')) = 1;
distance = letterindex.*levelindex(2,:);
distance = distance(indtemp+1);  % all monosac's
end