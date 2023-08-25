function [matched_gpcomp,matched_gpmass,matched_scanind] = ...
    findms1match(precmass,scannum,prot,protind,PTM,...
    searchscannum,scoreoptions)
% FINDMS1MATCH: find candidates that match precursor ion mass of spectra
%
% Syntax:
% [matched_gpcomp,matched_gpmass,matched_scanind] = ...
%     findms1match(precmass,scannum,prot,protind,PTM,...
%     searchscannum,scoreoptions)
%
% Input:
% precmass: n x 1 double. The precursor mass to be searched.
% scannum: n x 1 double. The scan number of corresponding spectra.
% prot: n x 1 cell array of m x 1 cell array of strings. The peptide 
%     backbones of candidate glycopeptides. See ASSEMBLEGP for details.
% protind: n x 1 numerical array. The index number of each protein. In
%     output "matched_gpcomp" these values will be used as protein ID.
% PTM: structure. The PTMs that will appear on peptides. See
%     SCOREALLSPECTRA and THEOPTMFRAG for details.
% searchscannum: m x 1 double. The scan numbers to be searched. If the scan
%     number in this input is not present in "scannum", no search will be
%     performed.
% scoreoptions: 1 x n structure. The options for scoring. 3 fields are
%     used: "ms1tol", "ms1tolunit" and "assemblegpoptions". See
%     SCOREALLSPECTRA for detail.
%
% Output:
% matched_gpcomp: n x 1 cell array of 1 x m double. The composition of
%     matched candidate. It contains the following information:
%     [Protein #, Peptide #, PTM 1 #, PTM 2 #,...].
%     Note that protein number is taken from input "protind".
% matched_gpmass: n x 1 double. The theoretical monosiotopic mass of
%     matched glycopeptide candidate.
% matched_scanind: n x 1 double. The index number of the scans that are
%     matched by a candidate. For example, suppose x = matched_scanind(i),
%     then for this glycopeptide - spectrum match, the composition of this
%     glycopeptide is matched_gpcomp{i}, its theoretical mass is
%     matched_gpmass(i), the matched spectrum's scan number is
%     searchscannum(x).
%
% Note:
% N/A
%
% Example:
% N/A
%
% Children function: 
% ASSEMBLEGP  PERFORMMS1MATCH
%

sectionsize = 100000;  % data block size, will be explained later
ptmseq = PTM.ptmseq;
ptmmass = PTM.ptmmass;
ptmmassfix = PTM.ptmmassfix;
ptmuniind = PTM.ptmuniind;

agpopt.mode = scoreoptions.assemblegpmode;
agpopt.maxOglyonpep = scoreoptions.maxOglyonpep;
% Build glycopeptides by combining peptide backbones with PTMs.
[gpmass,gpcomp] = assemblegp(prot,ptmmass,ptmmassfix,ptmuniind,ptmseq,agpopt);
for ii = 1:length(gpcomp)
    gpcomp{ii}(1) = protind(gpcomp{ii}(1));
end
ms1tol = scoreoptions.ms1tol;
ms1tolunit = scoreoptions.ms1tolunit;

% Below, the program is designed to minimize memory use.
% Since the number of spectrum - candidate matches cannot be forseen, the
%     size of result's storage will change constantly. In MATLAB it causes
%     serious RAM waste. However, by reducing the number of size changes the
%     waste can be attenuated.
matched_gpcomp = {};  % Initialize final output
matched_gpmass = [];
matched_scanind = [];
currentind = 0;  % initialize write index
% Instead of adding one item each time, we build temporary variables which
%     can store a large amount of results. System will allocate RAM for
%     them so writing results to these variables is much more efficient.
% When temporary variables are full, the content is transferred to final
%     output. New, empty temporary variables will be built, ready for a new
%     cycle.
% In this way, variable size changes is minimized.
temp_matched_gpcomp = cell(sectionsize,1);
temp_matched_gpmass = zeros(sectionsize,1);
temp_scanind = zeros(sectionsize,1);
for ii = 1:length(searchscannum)
    thisprecmass = precmass(scannum == searchscannum(ii));
    if ~isempty(thisprecmass)
        switch upper(ms1tolunit)
            case 'PPM'
                goodcandidateind = abs((gpmass - thisprecmass)/thisprecmass) * 1e6 <= ms1tol;
            case 'DA'
                goodcandidateind = abs(gpmass - thisprecmass) <= ms1tol;
        end
        if any(goodcandidateind)
            writegpcomp = gpcomp(goodcandidateind);
            writegpmass = gpmass(goodcandidateind);
            for jj = 1:length(writegpcomp)
                currentind = currentind + 1;
                if currentind > sectionsize  % when temp. var. is full
                    % transfer data to output
                    matched_gpcomp = [matched_gpcomp;temp_matched_gpcomp];
                    matched_gpmass = [matched_gpmass;temp_matched_gpmass];
                    matched_scanind = [matched_scanind;temp_scanind];
                    % build fresh, empty temp. var.
                    temp_matched_gpcomp = cell(sectionsize,1);
                    temp_matched_gpmass = zeros(sectionsize,1);
                    temp_scanind = zeros(sectionsize,1);
                    % reset write index
                    currentind = 1;
                end
                temp_matched_gpcomp{currentind} = writegpcomp{jj};
                temp_matched_gpmass(currentind) = writegpmass(jj);
                temp_scanind(currentind) = ii;
            end
        end
    end
end
% when calculation is finished, some data may be left in temp. var.
% they should be transfered to output.
temp_matched_gpcomp = temp_matched_gpcomp(temp_scanind > 0);
temp_matched_gpmass = temp_matched_gpmass(temp_scanind > 0);
temp_scanind = temp_scanind(temp_scanind > 0);

matched_gpcomp = [matched_gpcomp;temp_matched_gpcomp];
matched_gpmass = [matched_gpmass;temp_matched_gpmass];
matched_scanind = [matched_scanind;temp_scanind];
end