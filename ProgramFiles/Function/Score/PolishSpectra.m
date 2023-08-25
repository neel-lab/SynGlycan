function CleanSpectra=PolishSpectra(Spectra, CutOffMed, FracMax)
%POLISHSPECTRA: Remove small peaks corresponding to instrument noise from 
% input 'Spectra'. These noisy peaks have intensity < 'CutOffMed'*median
% provided they have intensity <'FracMac'*max Intensity (By default 
% CutOffMed=2 and FracMax=0.02)
%
% Syntax: 
% CleanSpectra=PolishSpectra(Spectra, CutOffMed, FracMax)
%
% Input:
% Spectra: n x 2 numerical array, the original spectrum;
% CutOffMed: below how many times of median intensity should the peak be
% removed.
% FracMax: below what fraction of the maximum intensity should the peak be
% removed.
%
% Output:
% CleanSpectra: spectrum with noise removed
% 
% Note:
% Only those peaks whose intensity below both criteria will be removed.
%
% Children function:
% N/A
% 
% See also:
% removePrecursorIon 

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 11/10/2020

IntMedian=median(Spectra(:,2));  % Median of all peak intensities
IntMax=max(Spectra(:,2));  % Base Peak
threshold = min(CutOffMed*IntMedian,FracMax*IntMax);  % Noise filtering threshold is the lesser of 2 thresholds
CleanSpectra_int = Spectra(Spectra(:,2) >= threshold,:);  % Remove the peaks below threshold
CleanSpectra=sortrows(CleanSpectra_int,-2); % CleanSpectra is rank ordered 
% based on intensity since this helps output
end