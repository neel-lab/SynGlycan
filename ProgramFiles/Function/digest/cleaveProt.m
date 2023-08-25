function Prod = cleaveProt(Prot,Enzyme,MissedMax,MinPepLen,MaxPepLen)
% CLEAVEPROT: digest protein with proteases, ignore effect of modifications.
%
% Syntax:
% Prod = cleaveProt(Prot,Enzyme,MissedMax,MinPepLen,MaxPepLen)
%
% Input:
% Prot: 1 x n string. The amino acid sequence of a protein. Modified amino
%     acids are replaced by single letter codes.
% Enzyme: 1 x n cell array of strings. Proteases' names.
% MissedMax: double. Maximum number of missed cleavages allowed.
% MinPepLen: double. Minimum length of digested peptide allowed.
% MaxPepLen: double. Maximum length of digested peptide allowed.
%
% Output:
% Prod: structure. Digested peptides and their details. Fields are:
%     pep: 1 x n cell array of strings. Digested peptides with modified
%         amino acids replaced by single letter codes.
%     start: 1 x n double. Starting position of peptide in protein.
%     fin: 1 x n double. Ending position of peptide in protein.
%     short: (obsolete) 1 x n cell array of strings. Digested peptides with
%         modifications.
%     start_nomed: (obsolete) 1 x n double. Starting position of peptide in
%         protein without modification.
%     fin_nomod: 1 x n double. Ending position of peptide in
%         protein without modification.
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
% See Also:
% DIGESTSGP

% Author: Sriram Neelamegham

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 04/01/2020

% Check if input of enzyme can be found in GlycoPAT local protease database
proteasedb = Protease.mklocaldb;
if iscell(Enzyme)
    for k = 1:length(Enzyme)
        if ~proteasedb.isKey(upper(Enzyme{k}))
            error('MATLAB:GlycoPAT:UNSUPPORTEDENZYME','UNSUPPORTED ENZYME TYPE');
        end
    end
else
    if ~proteasedb.isKey(upper(Enzyme))
        error('MATLAB:GlycoPAT:UNSUPPORTEDENZYME','UNSUPPORTED ENZYME TYPE');
    end
end

[protmodstart,protmodend] = regexp(Prot,'<\w>','start','end');
protmodlen = protmodend - protmodstart + 1;

CleavagePos = [];
% Retrieve protease information, then search for cleavage sites
if iscell(Enzyme)
    for k = 1:length(Enzyme)
        if(proteasedb.isKey(upper(Enzyme{k})))
            ezncleaveexpr = proteasedb(upper(Enzyme{k}));
            
            % Because regular expression in protease database recognizes
            %     only upper case letters and modified AA are lowercase,
            %     cleavage at modified AA is disabled
            CleavagePos = [CleavagePos,regexp(Prot, ezncleaveexpr)];
        else
            error('MATLAB:GlycoPAT:UNSUPPORTEDENZYME','UNSUPPORTED ENZYME TYPE');
        end
    end
else
    if proteasedb.isKey(upper(Enzyme))
        ezncleaveexpr =  proteasedb(upper(Enzyme));
        CleavagePos  = [CleavagePos,regexp(Prot, ezncleaveexpr)];
    else
        error('MATLAB:GlycoPAT:UNSUPPORTEDENZYME','UNSUPPORTED ENZYME TYPE');
    end
end
CleavagePos  = [0 CleavagePos length(Prot)];  % Adds the first and last aacidPoss in Prot
CleavagePos  = unique(CleavagePos);  % Removes duplicate in case first and last are also digestion points

% Perform cleavage
Prod = [];
PepCount = 1;  % Digestion of peptide with given minimum length and allowed missed cleavage
for i = 0:MissedMax
    cleaveposlength = length(CleavagePos);
    for j = 1:(cleaveposlength-i-1)
        start = CleavagePos(j)+1;  % Peptide starting point: cleavage site + 1;
        finish = CleavagePos(j+i+1);  % Peptide ending point: next (or further, considering missed cleavage) cleavage site
        peptide_long = Prot(start:finish);
        modstart = [strfind(peptide_long,'<'),length(peptide_long)+1];
        modend = [0,strfind(peptide_long,'>')];
        peptide_short = '';
        for k = 1:length(modstart)
            peptide_short = [peptide_short,peptide_long(modend(k)+1:modstart(k)-1)];
        end
        if ((length(peptide_short) >= MinPepLen) && (length(peptide_short) <= MaxPepLen))
            Prod.pep(PepCount) = {Prot(start:finish)};  % convert digested peptide into cell
            Prod.short(PepCount) = {peptide_short};
            Prod.start(PepCount) = start;
            modlenfix = sum(protmodlen(start > protmodend));  % fix AA position shift caused by mod.
            Prod.start_nomod(PepCount) = start - modlenfix;
            Prod.fin(PepCount) = finish;
            Prod.fin_nomod(PepCount) = start - modlenfix + length(peptide_short) - 1;
            PepCount = PepCount + 1;
        end
    end
end
end