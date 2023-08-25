function varargout = glycanFrag(varargin)
%GLYCANFRAG: Break input glycan/glycopeptide at break points and return
%    fragments
%
% Syntax: 
% [NewGlyPep,Bion]=glycanFrag(SmallGlyPep,breakPt)
% [NewGlyPep,Bion,NGPmsind,Bionmsind]=glycanFrag(SmallGlyPep,breakPt,...
%     monosacind)
%
% Input: 
% SmallGlyPep: SGP format glycan/glycopeptide sequence
% breakPt: 1-D numerical array, the position where fragmentation should
%     occur.
% monosacind: n x 3 matrix, monosac. index and starting/ending position of
%     curly brackets. Each row refers to one monosac. 
%
% Output: 
% NewGlyPep: string, Y side of glycopeptide fragments.
% Bion: cell array of strings, B side of glycan fragments.
% NGPmsind: monosac. index of NewGlyPep, format is the same as input
%     "monosacind", contents correspond to monosac. that remains in
%     NewGlyPep.
% Bionmsind: cell array, monosac. index of Bion, same format as NGPmsind.
%
% Note:
% N/A
% 
% Example:
% (Below all example use the 3 input mode, result "NewGlyPep" and "glyBion"
%     will be the same as 2 input
% 1. for glycopeptide
% SmallGlyPep='GYLN{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}CT{n{h{s}}{n{h{s}{f}}}}R';
% breakPt=[9,13,52];
% [st,ed] = regexp(SmallGlyPep,'(?<={)[^{}]+[{}]*?','start','end');
% monosacind = [1:length(st);st;ed]';
% [NewGlyPep,Bion,NGPmsind,Bionmsind]=glycanFrag(SmallGlyPep,breakPt,...
%     monosacind)
% 
% NewGlyPep =
% 
%     'GYLN{n{n}}CT{n{h{s}}{n{h{f}}}}R'
% 
% 
% Bion =
% 
%   1×3 cell array
% 
%     {'{s}'}    {'{h{h}}'}    {'{h{h}{h{h{h}}{h{h}}}}'}
% 
% 
% NGPmsind =
% 
%      1     6     6
%      2     8     8
%     12    41    41
%     13    43    43
%     14    45    45
%     15    49    49
%     16    51    51
%     18    56    56
% 
% 
% Bionmsind =
% 
%   3×1 cell array
% 
%     {1×3 double}
%     {2×3 double}
%     {7×3 double}
%
% 2. for glycopeptide (with glycan and another modification on the
% same amino acid)
% SmallGlyPep='GYLN{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}CT<s>{n{h{s}}{n{h{s}{f}}}}R';
% breakPt=[9,13,51];
% [st,ed] = regexp(SmallGlyPep,'(?<={)[^{}]+[{}]*?','start','end');
% monosacind = [1:length(st);st;ed]';
% [NewGlyPep,Bion,NGPmsind,Bionmsind]=glycanFrag(SmallGlyPep,breakPt,...
% monosacind)
% 
% NewGlyPep =
% 
%     'GYLN{n{n}}CT<s>{n{h{s}}}R'
% 
% 
% Bion =
% 
%   1×3 cell array
% 
%     {'{n{h{s}{f}}}'}    {'{h{h}}'}    {'{h{h}{h{h{h}}{h{h}}}}'}
% 
% 
% NGPmsind =
% 
%      1     6     6
%      2     8     8
%     12    44    44
%     13    46    46
%     14    48    48
% 
% 
% Bionmsind =
% 
%   3×1 cell array
% 
%     {4×3 double}
%     {2×3 double}
%     {7×3 double}
%
% 3. for glycan
% SmallGlyPep='{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}';
% breakPt=[7,26];
% [st,ed] = regexp(SmallGlyPep,'(?<={)[^{}]+[{}]*?','start','end');
% monosacind = [1:length(st);st;ed]';
% [NewGlyPep,Bion,NGPmsind,Bionmsind]=glycanFrag(SmallGlyPep,breakPt,...
% monosacind)
% 
% NewGlyPep =
% 
%     '{n{n{h{h{h{h}}{h}}}}}'
% 
% 
% Bion =
% 
%   1×2 cell array
% 
%     {'{h}'}    {'{h{h{h}}}'}
% 
% 
% NGPmsind =
% 
%      1     2     2
%      2     4     4
%      3     6     6
%      7    17    17
%      8    19    19
%      9    21    21
%     10    25    25
% 
% 
% Bionmsind =
% 
%   2×1 cell array
% 
%     {1×3 double}
%     {3×3 double}
%
%See also UQFRAGION, JOINGLYPEP, COMPILEFRAGS,GLYCANFRAG, BREAKGLYPEP,
%MULTISGPFRAG. 

% Author: Sriram Neelamegham
% Date Lastly Updated: 8/11/14

if nargin == 2
    SmallGlyPep = varargin{1};
    breakPt = varargin{2};
    calcindex = false;
elseif nargin == 3
    SmallGlyPep = varargin{1};
    breakPt = varargin{2};
    monosacind = varargin{3}(:,1);
    monosacstart = varargin{3}(:,2);
    monosacend = varargin{3}(:,3);
    NGPmsind = [];
    Bionmsind = {};
    calcindex = true;
    tempmonosacind = monosacind;
end
nFrag=length(breakPt);
NewGlyPep=SmallGlyPep;
Bion=[];
for j=nFrag:-1:1   % start with outermost bond
    openBrac=1;
    closeBrac=0;
    glyBion='{';
    NewGlyPep(breakPt(j))='';
    bracketcount=0;
    while(openBrac~=closeBrac)
        if (NewGlyPep(breakPt(j))=='{')
            openBrac=openBrac+1;
            bracketcount = bracketcount + 1;
        end
        if (NewGlyPep(breakPt(j))=='}')
            closeBrac=closeBrac+1;
        end
        glyBion=[glyBion,NewGlyPep(breakPt(j))];
        NewGlyPep(breakPt(j))='';
    end
    Bion=[Bion,cellstr(glyBion)];
    if calcindex
        tempBionind = find(monosacstart > breakPt(j),1,'first');
        tempsgpind = tempBionind:tempBionind + bracketcount;
        tempBionind = tempmonosacind(tempsgpind);
        tempmonosacind(tempsgpind) = [];
        if any(tempBionind)
            Bionmsind = [Bionmsind;varargin{3}(tempBionind,:)];
        end
    end
end
if calcindex
    tempBionmsind = [];
    for i = 1:length(Bionmsind)
        tempBionmsind = [tempBionmsind;Bionmsind{i}(:,1)];
    end
    [~,diffind] = setdiff(monosacind,tempBionmsind);
   NGPmsind = varargin{3}(diffind,:);
end
if calcindex
    varargout{1} = NewGlyPep;
    varargout{2} = Bion;
    varargout{3} = NGPmsind;
    varargout{4} = Bionmsind;
else
    varargout{1} = NewGlyPep;
    varargout{2} = Bion;
end
end