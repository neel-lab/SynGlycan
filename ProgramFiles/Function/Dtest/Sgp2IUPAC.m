function IUPACstr = Sgp2IUPAC(SGP2_0)
%Sgp2IUPAC: Conver sgp 2.0 to IUPAC format.
%
%  Syntax:
%     IUPACstr = Sgp2IUPAC(SGP2_0)
%
%  Input:
%     SGP2_0: glycan in sgp 2.0 format
%  
%  Output:
%     IUPACstr: glycan in IUPAC format
%
%  Example: 
%     M5 = '{GlcNAc(b1-?){GlcNAc(b1-4){Man(b1-4)[{Man(a1-6)[{Man(a1-3)}]{Man(a1-6)}}]{Man(a1-3)}}}}';
%     IUPACstr = Sgp2IUPAC(M5);

% Author: Yusen Zhou
% Date Lastly Updated: 05/17/16

IUPACstr = '';
SGP      = regexprep(SGP2_0,'[\{,\}]','');

% Take out the unit block
UnitPat   = '(\((.*?)\)[A-Z_a-z_0-9]*)|(\[)|(\])';
UnitBlock = regexp(SGP,UnitPat,'Match');
startIdx  = length(UnitBlock);

% Reorganize in IUPAC format
for i = 1 : length(UnitBlock)
    ithBlock = UnitBlock{startIdx};
    if(ithBlock==']')
        IUPACstr = [IUPACstr '['];
    elseif(ithBlock=='[')
        IUPACstr = [IUPACstr ']'];
    else
        LinkExp = '\((.*?)\)';
        LinkInfo = regexp(ithBlock,LinkExp,'Match');
        ithBlock = regexprep(ithBlock,LinkExp,'');
        ithBlock = [ithBlock LinkInfo{1}];
        IUPACstr = [IUPACstr ithBlock];
    end
    startIdx = startIdx-1;
end
end