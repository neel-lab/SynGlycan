function mass =  AAMW(whichAA)
% AAMW: a quick way to get amino acid residue weight, without the OH and H.
% 
% Syntax:
% mass =  AAMW(whichAA)
% 
% Input:
% whichAA: string, single letter code of amino acid;
% 
% Output:
% mass: double, mass of amino acid residue
% 
% Note:
% 1 amino acid each time
%
% Example:
% AAlist = {'G','A','V','L','I','M','F','W','P','S','T','C','Y','N','Q','D','E','K','R','H'};
% for ii = 1:length(AAlist)
%     mass = AAMW(AAlist{ii});
%     disp(AAlist{ii});
%     disp(mass);
% end
%
% (Result is not shown here)
%
% Children function: 
% N/A
% 

switch whichAA
    case 'G'
        mass = 57.0214635; %Gly
    case 'A'
        mass = 71.0371136; %Ala
    case 'S'
        mass = 87.0320282; %Ser
    case 'P'
        mass = 97.0527637; %Pro
    case 'V'
        mass = 99.0684137; %Val
    case 'T'
        mass = 101.047678; %Thr
    case 'C'
        mass = 103.009184; %Cys
    case 'L'
        mass = 113.084064; %Leu
    case 'I'
        mass = 113.084064; %Ile
    case 'N'
        mass = 114.042927; %Asn
    case 'D'
        mass = 115.0269429; %Asp
    case 'Q'
        mass = 128.0585771; %Gln
    case 'K'
        mass = 128.0949627; %Lys
    case 'E'
        mass = 129.0425929; %Glu
    case 'M'
        mass = 131.0404844; %Met
    case 'H'
        mass = 137.0589113; %His
    case 'F'
        mass = 147.0684137; %Phe
    case 'R'
        mass = 156.1011103; %Arg
    case 'Y'
        mass = 163.0633284; %Tyr
    case 'W'
        mass = 186.0793126; %Trp
end
end