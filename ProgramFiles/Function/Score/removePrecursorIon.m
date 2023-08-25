function [SpectraB,PrecIonInt] = removePrecursorIon(SpectraA,MH,z,MS2tol,MS2tolUnit)
% REMOVEPRECURSORION: Remove Precursor Ion Peaks from MS2 Spectrum 
%
% Syntax: 
% [SpectraB,PrecIonInt] = removePrecursorIon(SpectraA,MH,z,MS2tol,MS2tolUnit) 
%
% Input:
% SpectraA: MS2 Spectrum
% MH: M+H of precursor ion
% z: charge state of precursor ion
% MS2tol: MS2 tolerance
% MS2tolUnit: MS2 tolerance unit, either 'Da' or 'ppm'
% 
% Output:
% SpectraB: SpectraA without precursor ion peak 
% PrecIonInt: Total intensity of precursor ion isotopic peaks, including all charge states.
% 
% Example:
% N/A
% 
% See also:
% POLISHSPECTRA  SPECTRADENOISING

SpectraB = SpectraA;  % Copy of input, this variable is where operations take place
Hmass = 1.007825032;  % Mass of H+
PrecIonInt = 0;
maxisotopes = 5;  % Consider 5 isotopic peaks with higher m/z than monoisotopic 
%     peaks. This number does not include the monoisotopic form.
%     Different charge states were considered as well. 

% Search all possible charge states of precursor ion (1 through reported charge state of precursor ion)
for ii = 1: z
    zz = ii;  % Charge state to search
    mzlist = [];  % Receiver of identified precursor ion peaks
    for jj = 0:maxisotopes
        pMass = (MH+jj*Hmass+(zz-1)*Hmass)/zz;
        if strcmpi(MS2tolUnit,'Da')  % If MS2 tolerence's unit is Da
            tempmzlist=find(abs(SpectraB(:,1)-pMass)<=MS2tol);
        elseif(strcmpi(MS2tolUnit,'ppm'))  % If MS2 tolerence's unit is ppm
            tempmzlist=find(abs(SpectraB(:,1)-pMass)<=MS2tol/1e6*pMass);
        end
        mzlist = [mzlist;tempmzlist];  % Indices of peaks identified as precursor ion
    end
    precpks = SpectraB(mzlist,:);
    PrecIonInt = PrecIonInt + sum(precpks(:,2));  % Output: precursor ion total intensity
    SpectraB(mzlist,:)=[];  % Output: remove precursor ions
end

end