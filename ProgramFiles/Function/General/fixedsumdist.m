function array = fixedsumdist(tgt,mins,maxs)
% FIXEDSUMDIST: find all combinations that sum up to a given value, e.g.
% given 3 digits, [0,3], get all combinations that sums up to 5

%% STUPID WAY TO DO THE JOB, BETTER BE FIXED
% - SHOULD USE DYNAMIC PROGRAMMING
if length(mins) ~= length(maxs)
    error('CHECK INPUT');
end
mins = reshape(mins,1,[]);
maxs = reshape(maxs,1,[]);
replacenums = zeros(size(mins));
for i = 1:length(replacenums)
    replacenums(i) = maxs(i) - mins(i) + 1;
end
allcombi = getallcombi(replacenums);
temparray = allcombi - 1 + repmat(mins,size(allcombi,1),1);
array = temparray(sum(temparray,2) == tgt,:);

end