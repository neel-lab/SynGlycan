function combinations = combigroups(varargin)
% COMBIGROUPS: quickly get all combination of groups of number
% 
% Syntax:
% combinations = combigroups(grp1, grp2, grp3,...)
%
% Input:
% grp: m x n numerical array. Each row is a subgroup. After combining, the
%     sequence inside the subgroup remains.
%
% Output:
% combinations: all possible combinations of subgroups.
%
% Note:
% N/A
%
% Example:
% grp1 = [0 1;1 0];
% grp2 = 3;
% grp3 = magic(3)
% combinations = combigroups(grp1, grp2, grp3)
% 
% grp3 =
% 
%      8     1     6
%      3     5     7
%      4     9     2
% 
% 
% combinations =
% 
%      0     1     3     8     1     6
%      0     1     3     3     5     7
%      0     1     3     4     9     2
%      1     0     3     8     1     6
%      1     0     3     3     5     7
%      1     0     3     4     9     2
%
% Children function: 
% N/A
%

if isnumeric(varargin{1})
    numofinputs = nargin;
    inputs = varargin;
elseif iscell(varargin{1})
    inputs = varargin{1};
    numofinputs = length(inputs);
end
numofelements = zeros(1,numofinputs);
outputwidth = 0;
for i = 1:numofinputs
    numofelements(i) = size(inputs{i},1);
    outputwidth = outputwidth + size(inputs{i},2);
end
combipattern = getallcombi(numofelements);
combinations = zeros(size(combipattern,1),outputwidth);
for i = 1:size(combipattern,1)
    tempcombi = [];
    for j = 1:size(combipattern,2)
        tempcombi = [tempcombi,inputs{j}(combipattern(i,j),:)];
    end
    combinations(i,:) = tempcombi;
end

end