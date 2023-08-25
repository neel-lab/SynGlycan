function structfeatures = getglycanstructfeatures(glyseq)
%
% Classifies N-glycan into 1 of 17 types based in list below
% 
% Input: sgp1.0 string called glyseq
% Output: 17X1 numerical array. Number in each location represents a
% structural feature based in the list that follows
%
% Example:
% feature= getglycanstructfeatures('{n{f}{n{h{h{n{h}}{n{f}{h{s}}}}{h{h}{h}}}}}')
glyseq = char(glyseq);

is_high_mannose = 0;  % 1
is_pauci_mannose = 0;  % 2
is_hybrid = 0;  % 3
is_complex = 0;  % 4
is_fucosylated = 0;  % 5
not_fucosylated = 0; % 5
is_core_fucosylated_only = 0;  % 6
is_terminal_fucosylated_only = 0;  % 7
is_both_core_and_terminal_fucosylated = 0;  %8
is_terminated_by_Neu5Ac = 0;  % 9
is_terminated_by_Gal = 0;  % 10
is_terminated_by_GlcNAc = 0;  % 11
is_mono_antennary = 0;  % 12
is_bi_antennary = 0;  % 13
is_tri_antennary = 0;  % 14
is_tetra_antennary = 0;  % 15
is_other_antennary = 0;  % 16
is_bisecting = 0;  % 17
is_lewis = 0;
is_slewis = 0;
is_lacdiNAc = 0;

numh = length(strfind(glyseq,'h'));
numn = length(strfind(glyseq,'n'));
numf = length(strfind(glyseq,'f'));
nums = length(strfind(glyseq,'s'));
distance = getdistance(glyseq);
monosac = regexp(glyseq,'[hnsf]','match');
hdist = distance(strcmpi(monosac,'h'));
ndist = distance(strcmpi(monosac,'n'));
fdist = distance(strcmpi(monosac,'f'));
bondmap = zeros(length(distance));
readind = 1;

while readind < length(distance)
    if (distance(readind + 1) - distance(readind)) == 1
        bondmap(readind,readind + 1) = 1;  % consecutive numbers indicate bond
        readind = readind + 1;
    elseif (distance(readind + 1) - distance(readind)) < 1  % if chain is broken, go back to find its fork point
        thisind = distance(readind + 1);  % where it's broken
        itsforkpt = find(distance(1:readind) == thisind - 1,1,'last');  % where is the fork point
        bondmap(itsforkpt,readind + 1) = 1;  % mark this bond
        readind = readind + 1;  % keep going on
    end
end

if numn == 2 && numh > 4
    is_high_mannose = 1;  % high-mannose
elseif numn == 2 && numh > 0
    is_pauci_mannose = 1;  % pauci-mannose
else
    if sum(hdist == 5) > 0
        is_hybrid = 1;  % hybrid
    else
        if any(ndist ==5)
            is_complex = 1;  % complex
            % only complex glycans have antennae
            numantennae = sum(ndist ==5);
            switch numantennae
                case 1
                    is_mono_antennary = 1;
                    nant=1;
                case 2
                    is_bi_antennary = 1;
                    nant=2;
                case 3
                    is_tri_antennary = 1;
                    nant=3;
                case 4
                    is_tetra_antennary = 1;
                    nant=4;
                otherwise
                    is_other_antennary = 1;
                    nant=1;
            end
        end
    end
end
% fucosylated
if any(numf)
    is_fucosylated = 1;  % has fucose
    if any(fdist == 2)
        if any(fdist > 2)  % core fucosylated + terminal
            is_both_core_and_terminal_fucosylated = 1;
        else  % core fucosylated only
            is_core_fucosylated_only = 1;
        end
    else  % terminal fucosylated only
        is_terminal_fucosylated_only = 1;
    end
else
    not_fucosylated = 1;
end
if is_complex > 0
    allnind = find(strcmpi(monosac,'n'));
    antennaeroot = allnind(ndist == 5);
    for i = 1:length(antennaeroot)
        childrenind = glytreetracker(bondmap,antennaeroot(i),[],'down');
        childrenmonosac = monosac(childrenind);
        % Here the logic is: if a branch has Neu5Ac, it's terminated by Neu5Ac
        % otherwise if a branch has Gal, it's terminated by Gal
        % otherwise if a branch has GlcNAc, it's terminated by GlcNAc
        if ismember('s',childrenmonosac)
            is_terminated_by_Neu5Ac = is_terminated_by_Neu5Ac + 1;
        elseif ismember('h',childrenmonosac)
            is_terminated_by_Gal = is_terminated_by_Gal + 1;
        elseif ismember('n',childrenmonosac)
            is_terminated_by_GlcNAc = is_terminated_by_GlcNAc + 1;
        end
    end
    is_terminated_by_Neu5Ac=is_terminated_by_Neu5Ac/length(antennaeroot);
    is_terminated_by_Gal=is_terminated_by_Gal/length(antennaeroot);
    is_terminated_by_GlcNAc=is_terminated_by_GlcNAc/length(antennaeroot);
end
% bisecting
if any(ndist == 4)
    is_bisecting = 1;
end
% lewis, sialyl-lewis,lacdiNAc
if any(ndist == 6)
    is_lacdiNAc = 1;
end
if any(strfind(glyseq,'{n{f}{h}}')) || any(strfind(glyseq,'{n{h}{f}}'))
    is_lewis = 1;
end
if any(strfind(glyseq,'{n{f}{h{s}}}')) || any(strfind(glyseq,'{n{h{s}}{f}}'))
    is_slewis = 1;
end
structfeatures = [is_high_mannose,is_pauci_mannose,is_hybrid,...
    is_complex,not_fucosylated,is_core_fucosylated_only,...
    is_terminal_fucosylated_only,is_both_core_and_terminal_fucosylated,is_terminated_by_Neu5Ac,...
    is_terminated_by_Gal,is_terminated_by_GlcNAc,is_mono_antennary,...
    is_bi_antennary,is_tri_antennary,is_tetra_antennary,...
    is_other_antennary,is_bisecting,is_lacdiNAc,is_lewis,is_slewis];
end