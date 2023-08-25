function Top10 = calculatetop10match(spectrum,peakismatched,windowidth,Da_min)
% CALCULATETOP10MATCH: calculate the Top 10 value of spectrum-candidate
% pair
%
% Syntax:
% Top10 = calculatetop10match(spectrum,peakismatched,windowidth,Da_min)
%
% Input:
% spectrum: experimental spectrum. spectrum has been sorted by descending
% intensity
% peakismatched: n x 1 cell or logical array, whether a fragment peak is
%     matched (and which peak, if input is a cell array).
% windowidth: when a matched peak is among the top 10, other peaks inside
%     this window, if matched, will not be counted as top 10. By default we
%     exclude peaks within +/1 1Da s that isotope peaks are not matched.
% Da_min: only peaks above this m/z threshold will be considered in
%     calculation. By default Da_min is set to 400 Da in order to exlcude B-ions in HCD
% 
% Output:
% Top10: Number of the top10 that are matched (numerical value)
%
% Note:
% N/A
%
% Example:
% N/A
%
% Children function: 
% N/A
%
% 
if iscell(peakismatched)
    peakismatched = ~cellfun(@isempty,peakismatched);
end

halfwidth = windowidth/2;
groupnum = zeros(size(spectrum,1),1);
currentgrpnum = 1;
groupnum(1) = currentgrpnum;
if ~isempty(spectrum)
    % pre-treatment: based on peak intensity and user set half width, put
    % the peaks in groups. Top 10 calculation is performed onto groups,
    % rather than individual peaks
    % this calc has nothing to do with matching results.
    boundaries = [spectrum(:,1) - halfwidth,spectrum(:,1) + halfwidth];
    for i = 2:size(spectrum,1)
        % whether this peak is inside a group defined previously
        isinrange = spectrum(i,1) >= boundaries(1:i-1,1) & spectrum(i,1) <= boundaries(1:i-1,2);
        if any(isinrange)  % find out which group, assign this peak with that group number
            tempgroupnum = groupnum(1:i);
            matchedgrpnum = tempgroupnum(isinrange);
            groupnum(i) = matchedgrpnum(1);
        else  % assign this peak with a new group number
            currentgrpnum = currentgrpnum + 1;
            groupnum(i) = currentgrpnum;
        end
    end
end
aboveDa_min = spectrum(:,1) >= Da_min;
peakismatched_above = peakismatched(aboveDa_min);
groupnum_above = groupnum(aboveDa_min);
grpcounter = 0;
index = 1;
% among the first 10 groups, how many of them contain at least 1 matched
% peak. this number is the TOP10 we need
while grpcounter < 10 && index <= length(groupnum_above)
    prevgrps = groupnum_above(1:index-1);
    if ~ismember(groupnum_above(index),prevgrps)
        grpcounter = grpcounter + 1;
    end
    index = index + 1;
end
% if a peak is low but in same group with a high peak, if high peak is not
% matched, even if the low peak is matched, it will not be counted as one
% of TOP10. This is to prevent overestimation

temp_groupnum_above = groupnum_above(1:index-1);
temp_peakismatched_above = peakismatched_above(1:index-1);
matchedgrps = temp_groupnum_above(temp_peakismatched_above);
Top10 = length(unique(matchedgrps));

end