function IUPACout = SGP2toIUPAC(SGP2in)
indtemp = strfind(SGP2in,'{');
levelindex = zeros(2,length(SGP2in));
indtemp2 = strfind(SGP2in,'}');
if ~ismember(1,indtemp)
    error('Check input sequence, DrawGlycan cannot read this.');
end
levelindex(1,indtemp) = 1;
levelindex(1,indtemp2) = -1;
levelindex(2,1) = 1;
for ii = 2:size(levelindex,2)
    levelindex(2,ii) = levelindex(2,ii-1) + levelindex(1,ii);
end
[monosac,monosacstart] = regexp(SGP2in,'\(.*?\)[A-z0-9]+','match','start');
for ii = 1:length(monosac)
    temp = monosac{ii};
    [linkstart,linkend] = regexp(temp,'\(.*?\)','start','end');
    monosac{ii} = [temp(linkend + 1:end),temp(linkstart:linkend)];
end
IUPACout = '';
%         insidebranch = false;
bracketlvlstack = [];
for ii = length(monosac):-1:1
    inhibitor = false;
    detectorind = monosacstart(ii) - 2;
    if detectorind > 0
        if strcmpi(SGP2in(detectorind),'}')  % incoming branching point, insert square bracket
            bracketlvlstack = [bracketlvlstack;levelindex(2,detectorind)];
            IUPACout = [IUPACout,monosac{ii},'['];
            inhibitor = true;
        else  % regular
            IUPACout = [IUPACout,monosac{ii}];
        end
        if ~isempty(bracketlvlstack) && levelindex(2,detectorind) == bracketlvlstack(end) && ~inhibitor
            IUPACout = [IUPACout,']'];
            bracketlvlstack(end) = [];
        end
    else  % very first monosac
        IUPACout = [IUPACout,monosac{ii}];
    end
end
left = count(IUPACout,'[');
right = count(IUPACout,']');
if right ~= left
    while left ~= right
        if strcmpi(IUPACout(end-32),'M')
            IUPACseqleft = IUPACout(1:(end-33));
            IUPACseqright= IUPACout((end-32):end);
        else
            IUPACseqleft = IUPACout(1:(end-35));
            IUPACseqright= IUPACout((end-34):end);
        end
        IUPACout = strcat(IUPACseqleft,']',IUPACseqright);
        left = count(IUPACout,'[');
        right = count(IUPACout,']');
    end
end
end