function out = sgp1to2a(in,PTMid)
%% sgp1 to sgp2
if PTMid == 1
    [msstart] = regexp(in,'(?<=[{}])[a-z]+(?=[{}])','start');
    index = getmsdist(in);
    if ~contains(in,'{h{h')
        type='O';
    else
        type='N';
    end
    out = in;
    for i = length(msstart):-1:1
        thisms = in(msstart(i));
        switch thisms
            case 'n'
                newms = '(??-?)GlcNAc';  %N-Glycan
                if type=='O' && index(i) == 1
                    newms = '(??-?)GalNAc';
                end
%                 newms = '(??-?)GalNAc';  %O-Glycan
            case 'h'
                newms = '(??-?)Gal';
                if type=='N' && (index(i) <= 5 && index(i) >= 3)
                  newms = '(??-?)Man';  %N-Glycan
%                if index(i) <= 5 && index(i) >= 3
%                    newms = '(??-?)Man';  %N-Glycan
%                     newms = '(??-?)Gal';  %O-Glycan
%                else
%                    if (msstart(i) - 2)>0
%                    if strcmpi(in(msstart(i) - 2),'h')
%                        newms = '(??-?)Man';  %N-Glycan
%                       newms = '(??-?)Gal';  %O-Glycan
%                    else
%                        newms = '(??-?)Gal';
%                    end
%                    else
%                        newms = '(??-?)Gal';
%                    end
%}
                end
            case 's'
                newms = '(??-?)Neu5Ac';
            case 'f'
                newms = '(??-?)Fuc';
            case 'g'
                newms = '(??-?)Neu5Gc';
            case 'x'
                newms = '(??-?)Xyl';
            case 'l'
                newms = '(??-?)Deoxynonulosonate';
        end
        out = [out(1:msstart(i)-1),newms,out(msstart(i)+1:end)];
    end
    if isempty(msstart) && strcmpi(in(1),'<')
        out = ['{(??-?)',out(2:end-1),'}'];
    end
elseif PTMid == 0
    if strcmpi(in(1),'<') && strcmpi(in(end),'>')
        out = ['{(??-?)',in(2:end-1),'}'];
    end
end
end
function distance = getmsdist(thisgly)
indtemp = strfind(thisgly,'{');
levelindex = zeros(2,length(thisgly));
indtemp2 = strfind(thisgly,'}');
levelindex(1,indtemp) = 1;
levelindex(1,indtemp2) = -1;
levelindex(2,1) = 1;
for j = 2:size(levelindex,2)
    levelindex(2,j) = levelindex(2,j-1) + levelindex(1,j);
end
wholestr = zeros(1,length(thisgly));
wholestr(indtemp) = 1;
wholestr(indtemp2) = 1;
readind = 1;
allms = cell(sum(wholestr)/2,1);
writeind = 1;
tempstr = '';
while readind <= length(wholestr)
    if wholestr(readind) ~= 1  % gather character to form the monosac. string
        tempstr = [tempstr,thisgly(readind)];
    else
        if ~isempty(tempstr)
            bondinfo = strfind(tempstr,')');
            if ~isempty(bondinfo)  % using parenthesis
                allms{writeind} = tempstr(max(bondinfo) + 1:end);
            else  % No parenthesis means there's only monosac.
                allms{writeind} = tempstr;
            end
            tempstr = '';
            writeind = writeind + 1;
        end
    end
    readind = readind + 1;
end
%% calculate distance of each monosac
letterindex = zeros(1,length(thisgly));
letterindex(regexp(thisgly,'[^{}]')) = 1;
distance = letterindex.*levelindex(2,:);
distance = distance(indtemp+1);  % all monosac's
end