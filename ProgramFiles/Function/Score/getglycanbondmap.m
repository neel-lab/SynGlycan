function bondmap = getglycanbondmap(thisgly)
% GETGLYCANBONDMAP: get connection table for glycan structure
%
% Syntax:
% bondmap = getglycanbondmap(thisgly)
% 
% Input:
% thisgly: SGP 1.0 sequence of the glycan.
% 
% Output:
% bondmap: n x n numerical array describing glycan tree structure.
%     bondmap(m,n) = 1 means there is a linkage between mth and nth
%     monosac. 0 means no connection. Indices of monosac. is given by
%     counting from the left side of its SGP sequence. By default n > m. 
% 
% Note:
% bondmap is unidirectional, an upper triangular matrix.
%
% Example:
% bondmap = getglycanbondmap('{n{f}{h{s}}}')
% 
% bondmap =
% 
%      0     1     1     0
%      0     0     0     0
%      0     0     0     1
%      0     0     0     0
% 
% Children function: 
% N/A


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
end