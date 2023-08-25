function allcombi = getallcombi(maxnums)
% GETALLCOMBI: By assigning the maximum of each element, return all possible
% combinations.
%
% Syntax:
% allcombi = getallcombi(maxnums)
%
% Input:
% maxnums: 1 x n numerical array. The maximum value of each number in the
% output. 
%
% Output:
% allcombi: m x n numerical array. n equals to input width. Each element
% varies from 1 to the number at the corresponding position in the input.
%
% Note:
% All inputs must be non-negative integer.
%
% Example:
% allcombi = getallcombi([2 0 3 0 1])
% 
% allcombi =
% 
%      1     0     1     0     1
%      1     0     2     0     1
%      1     0     3     0     1
%      2     0     1     0     1
%      2     0     2     0     1
%      2     0     3     0     1
%
% Children function: 
% N/A
%

if ~isempty(maxnums)
    maxnums = reshape(maxnums,1,length(maxnums));
    inputiszero = maxnums == 0;
    maxnums(inputiszero) = 1;  % The above steps is to handle zero's, replace them with 1's then replace back
    totnum = prod(maxnums);
    allcombi = zeros(totnum,length(maxnums));
    dupfactor = 1;
    for i = length(maxnums):-1:1
        tempcolumn = 1:maxnums(i);
        tempcolumn = repmat(tempcolumn,dupfactor,1);
        tempcolumn = tempcolumn(:);
        allcombi(:,i) = repmat(tempcolumn,totnum/length(tempcolumn),1);
        dupfactor = dupfactor * maxnums(i);
    end
    allcombi(:,inputiszero) = 0;
else
    allcombi = [];
end
end