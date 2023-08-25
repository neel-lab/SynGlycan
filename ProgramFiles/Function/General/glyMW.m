function gmw = glyMW(varargin)
% GLYMW: Calculate monoisotopic MW of glycan (NOT M+H)
% 
% Syntax: 
%      glyFull = glyMW(glycanstring)
% 
%   Input:
%     glycanstring: Glycan string, optionally enclosed by {}
% 
%   Output: 
%     glyFull: monoisotopic glycan MW (NOT M+H)
%
% Examples: 
%  Example 1: 
%   for 'HexNAcHexNeuAc' (brackets are ignored)
%   >> glyMW('{n{h{s}}}')
%   >> 
%   Answer:
%       ans=674.2381756  
%
%  Example 2:
%   for NeuAc2,3Hex2,3[NeuAc2,6]HexNAc (linkage specificity mentioned in the text is ignored) 
%    >> glyMW('{n{s}{h{s}}}')
%   Answer:
%       ans = 965.333592
%
%  Example 3 (no brackets with parentheses):
%     HexNAcHex(NeuAc)
%   >> glyMW('nh(s)'); 
%    Answer:
%       ans = 674.2381756
%
%  Example 4:
%   for FucNeuAcHexHexNAc
%    >> slex='fshn'; 
%      glyMW(strcat(slex,'hn'))
%    Answer:
%       ans = 1185.4282804
%
%See also: SCOREGUI, BROWSEGUI, SCOREALLSPECTRA, GLYMZCALC, GLYFORMULA

% Author: Sriram Neelamegham, Gabrielle Chapman, and Gang Liu
% Date Lastly Updated: 10/27/22 by Gabrielle Chapman

% Note: Solution verification:
% Answers checked using GlycanMass

if (nargin>0)
    gly=varargin{1};
end

if (nargin>1)
    GlycanNew=varargin{2};   % This part of the code is not written but it will yield glycan molecular composition when done
    AnomericNew=varargin{3};
else
    GlycanNew = Glycan;
    AnomericNew = Anomeric;
end

format longg;
glymw=0;
if ischar(gly)
    [starting,ending]  = regexp(gly,'\d+(\.\d+)?');  % find number in string
    for k=1:length(starting)
        str=(gly(starting(k):ending(k)));
        glymw=glymw+str2double(str);
    end
    
    typesofmono = GlycanNew.glycanMSMap.keys;                     
    numtypesofmono = length(typesofmono);
    for i = 1 : numtypesofmono
        nummonoresidues = length(strfind(gly,typesofmono{i}));
        glymw = glymw+nummonoresidues*GlycanNew.glycanMSMap(typesofmono{i});
    end
    gmw=glymw+18.0105633;
else
    [starting,ending]  = regexp(gly.name,'\d+(\.\d+)?');  % find number in string
    for k=1:length(starting)
        str=(gly.name(starting(k):ending(k)));
        glymw=glymw+str2double(str);
    end
    form = gly.form;
    if any(strcmpi(form,{'Me','MY','MB','MI','MA','MBX','MXY'}))         % for all methylated ions
        typesofmono = GlycanNew.glycanMSMeMap.keys;
        numtypesofmono = length(typesofmono);
        for i = 1 : numtypesofmono
            nummonoresidues = length(strfind(gly.name,typesofmono{i}));
            glymw = glymw+nummonoresidues*GlycanNew.glycanMSMeMap(typesofmono{i});
        end
    elseif any(strcmpi(form,{'Ac','AY','AB','AI','AA','ABX','AXY'}))     % for all acetylated ions
        typesofmono = GlycanNew.glycanMSAcetylMap.keys;
        numtypesofmono = length(typesofmono);
        for i = 1 : numtypesofmono
            nummonoresidues = length(strfind(gly.name,typesofmono{i}));
            glymw = glymw+nummonoresidues*GlycanNew.glycanMSAcetylMap(typesofmono{i});
        end
    elseif any(strcmpi(form,{'OH','Y','B','I','A','BX','XY'}))      % for non-derivatized ions
        typesofmono = GlycanNew.glycanMSMap.keys;
        numtypesofmono = length(typesofmono);
        for i = 1 : numtypesofmono
            nummonoresidues = length(strfind(gly.name,typesofmono{i}));
            glymw = glymw+nummonoresidues*GlycanNew.glycanMSMap(typesofmono{i});
        end
    end                    
    redEnd=AnomericNew.mass(gly.ano);
    % non-reducing end addition
    if any(strcmpi(form,{'Me','MY','Ac','AY','OH','Y','XY','MXY','AXY'}))
        if any(strcmpi(form,{'MY','AY','Y','OH','XY','MXY','AXY'}))
            TotEnd=redEnd+1.0078246;    % Add +H at non-reducing end
        elseif any(strcmpi(form,{'Me'})) % Methyl
            TotEnd=redEnd+15.0234738;    % Add +CH3 at non-reducing end
        elseif any(strcmpi(form,{'Ac'})) % Acetyl
            TotEnd=redEnd+43.0183879;     % Add +COCH3 at non-reducing end
        end
    end
    if any(strcmpi(form,{'B','MB','AB','BX','MBX','ABX'}))
        if any(strcmpi(form,{'MB'}))  % Methyl
            TotEnd=15.0234738;        % Add +CH3 at non-reducing end
        elseif any(strcmpi(form,{'AB'}))  % Acetyl
            TotEnd=43.0183879;             % Add +COCH3 at non-reducing end
        else
            TotEnd=1.0078246;           % Add +H at non-reducing end
        end
        TotEnd=TotEnd-1.0078246;    % Add H at reducing end
    end
    if any(strcmpi(form,{'MI','AI','I','A','MA','AA'}))
        TotEnd=0; % loss of H at red-end & gain at other end compensate
    end
    gmw=glymw+TotEnd;
end
end