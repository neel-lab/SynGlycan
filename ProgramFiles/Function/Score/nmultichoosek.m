function combs = nmultichoosek(values, k)
% NMULTICHOOSEK: create combinations with repetitions, similar to
%     NCHOOSEK but with repeats.
% Example:
% combs = nmultichoosek([0,2,9],2)
% 
% combs =
% 
%      0     0
%      0     2
%      0     9
%      2     2
%      2     9
%      9     9
% 


% Quoted from stackoverflow.com/users/3139711/knedlsepp

if numel(values)==1 
    n = values;
    combs = nchoosek(n+k-1,k);
else
    n = numel(values);
    combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
    combs = reshape(values(combs),[],k);
end