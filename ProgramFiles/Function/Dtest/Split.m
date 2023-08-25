function [MonoStr,MonoIndex,BracStr,BracIndex] = Split(IUPAC)
% SPLIT: split IUPAC sequence input to retrieve monosaccharide string,
%     position, brackets and their positions.
%
% Syntax:
% [MonoStr,MonoIndex,BracStr,BracIndex] = Split(IUPAC)
%
% Input:
% IUPAC: string, IUPAC sequence input.
%
% Output:
% MonoStr: 1 x n cell array of strings. Monosaccharide name, including
%     content enclosed in parenthesis.
% MonoIndex: 1 x n numerical array. Position of the first letter of
%     monosaccharide name.
% BracStr: 1 x n cell array of strings. Bracket symbols, either "[" or "]".
% BracIndex: 1 x n numerical array. Position of brackets.
%
% Note:
% Leading/trailing whitespaces are included in monosaccharide name
%     (including parenthesis), but not their positions. For example "  Gal ()  "
%     will be retrieved as "  Gal ()  " but "MonoIndex" will be 3.
% This program uses a "delimiter based" approach. Monosaccharides are
%     divided by parentheses as well as brackets, therefore anything in between
%     are treated as monosaccharide name.
%
% Example:
% [MonoStr,MonoIndex,BracStr,BracIndex] = Split('  Aa($%&&*)Bb[Ff () G] C [Hh()[]Ii()Jj] Dd()Ee ')
%
% MonoStr =
%
%     '  Aa($%&&*)'    'Bb'    'Ff ()'    ' G'    ' C '    'Hh()'    'Ii()'    'Jj'    ' Dd()'    'Ee '
%
%
% MonoIndex =
%
%      3    12    15    21    24    27    33    37    41    45
%
%
% BracStr =
%
%     '['    ']'    '['    '['    ']'    ']'
%
%
% BracIndex =
%
%     14    22    26    31    32    39
%
% Children function:
% N/A
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2020, Research Foundation for State University of New York. All rights reserved
%

% Block user costomized info to avoid confusion when locating local options
[optcontstart,optcontend] = regexp(IUPAC,DrawGlycanPara.regexp_optionvalue,'start','end');
temp_IUPAC = IUPAC;
for i = 1:length(optcontstart)
    temp_IUPAC(optcontstart(i):optcontend(i)) = blanks(optcontend(i)-optcontstart(i)+1);
end
% Search parenthesis
[posparenpairstart,posparenpairend] = regexp(temp_IUPAC,'\(.*?\)','start','end');
% Search brackets
[BracStr,BracIndex] = regexp(temp_IUPAC,'[\]\[]','match','start');
alldelimiters = [posparenpairstart,posparenpairend,BracIndex];
% Anything that's not delimiters are monosaccharide names
allmono = setdiff(1:length(temp_IUPAC),alldelimiters);
% Extract monosac. names by searching continuous segments
temp_allmonoind = find((diff(allmono) == 1) == 0);
allmonoind = [[1,temp_allmonoind+1];[temp_allmonoind,length(allmono)]];  % positions for monosac. names
tempallmonopos = [allmono(allmonoind(1,:));allmono(allmonoind(2,:))];
allparenpairpos = [posparenpairstart;posparenpairend];
if isempty(allparenpairpos)
    allparenpairpos = [0;0];
end
ind = 1;
allmonopos = [];
% "tempallmonopos" is a mixture of info inside and outside parenthesis,
%     they need to be reconstituted before output. Compare their starting
%     and ending locations with those of parenthesis will do the job. The
%     rule is: parenthesis always appear at the end of a monosac. name
%     string.
while ind <= size(tempallmonopos,2)
    thismonopos = tempallmonopos(:,ind);
    parenthesismatch = (thismonopos(2) + 1 == allparenpairpos(1,:));
    if any(parenthesismatch)
        thismonopos(2) = allparenpairpos(2,parenthesismatch);
        allmonopos = [allmonopos,thismonopos];
        ind = find(tempallmonopos(1,:) > thismonopos(2),1);
    else
        allmonopos = [allmonopos,thismonopos];
        ind = ind + 1;
    end
end
MonoStr = cell(1,size(allmonopos,2));
% Below is how "Note" section's first item works (monosac. name starting
%     location is not affected by whitespaces).
whitespacefix = zeros(1,length(MonoStr));
for i = 1:length(MonoStr)
    MonoStr{i} = IUPAC(allmonopos(1,i):allmonopos(2,i));
    whitespacefix(i) = find(double(MonoStr{i}) - 32 ~= 0,1);
end
MonoIndex = allmonopos(1,:);
MonoIndex = MonoIndex + whitespacefix - 1;
end