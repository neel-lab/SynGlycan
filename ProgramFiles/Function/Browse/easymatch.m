function foundPeaks = easymatch(spectrum,sgp,result,theofrag)
% EASYMATCH: find the information of each matched peak during the browsing
%     process.
% 
% Syntax:
% foundPeaks = easymatch(spectrum,sgp,result,theofrag)
% 
% Input:
% spectrum: n x 2 double. Spectrum.
% sgp: string. Glycopeptide sequence.
% result: structure. Information about matched peaks. Field
%     "peakmatchindex" is used. It contains the theoretical fragments
%     corresponding to each mathced peak.
% theofrag: 1 x nn structure. All theoretical fragments of candidate
%     glycopeptide.
% 
% Output:
% foundPeaks: 1 x n structure. All theoretical fragments that are matched.
%     Fields are:
%         original: string. SGP of glycopeptide.
%         sgp: string. SGP of fragment.
%         nmFrag: double. Number of non-glycan PTM cleavage(s) on this fragment.
%         npFrag: double. Number of glycan PTM cleavage(s) on this fragment.
%         ngFrag: double. Number of peptide cleavage(s) on this fragment.
%         mz: double. m/z of fragment.
%         type: string. Ion type of fragment.
%         charge: double. Charge state of fragment
%         mzTheo: double. 
%         mzExpt: double. 
%         ppmError: double. 
%         DaError: double. 
%         Intensity: double. 
%         peakIndex: double. 
%         iontype
%         iontype2
%         unitindex
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
% 

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 11/10/2020

foundPeaks = [];
pkmatchind = result.peakmatchindex;
fcount = 1;
for i = 1:length(pkmatchind)
    temppmi = pkmatchind{i};
    if ~isempty(temppmi)
        for j = 1:size(temppmi,1)
            foundPeaks(fcount).original=sgp;
            foundPeaks(fcount).sgp=theofrag(temppmi(j,1)).sgp;
            foundPeaks(fcount).nmFrag=theofrag(temppmi(j,1)).nmFrag;
            foundPeaks(fcount).npFrag=theofrag(temppmi(j,1)).npFrag;
            foundPeaks(fcount).ngFrag=theofrag(temppmi(j,1)).ngFrag;
            foundPeaks(fcount).mz=theofrag(temppmi(j,1)).mz;
            foundPeaks(fcount).type=theofrag(temppmi(j,1)).type;
            foundPeaks(fcount).charge=temppmi(j,2);
            foundPeaks(fcount).mzTheo=(theofrag(temppmi(j,1)).mz - 1.007825032)/temppmi(j,2) + 1.007825032;
            foundPeaks(fcount).mzExpt=spectrum(i,1);
            foundPeaks(fcount).ppmError=(spectrum(i,1)-foundPeaks(fcount).mzTheo)/foundPeaks(fcount).mzTheo*1e6;
            foundPeaks(fcount).DaError = spectrum(i,1)-foundPeaks(fcount).mzTheo;
            foundPeaks(fcount).Intensity=spectrum(i,2);
            foundPeaks(fcount).peakIndex = i;
            foundPeaks(fcount).iontype = theofrag(temppmi(j,1)).type;
            foundPeaks(fcount).iontype2 = theofrag(temppmi(j,1)).type;
            foundPeaks(fcount).unitindex = theofrag(temppmi(j,1)).unitindex;
            fcount = fcount + 1;
        end
    end
end
end