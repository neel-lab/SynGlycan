function spectra = spectradenoising(spectra,datafragmodes,datacharges,...
    datamonomasses,ms2tols,ms2tolunits,denoisingoptions)
% SPECTRADENOISING: delete noise peaks in spectrum
%
% Syntax:
% spectra = spectradenoising(spectra,fragmodes,charges,...
%     monomasses,ms2tols,ms2tolunits,scoreoptions)
%
% Input:
% spectra: m x 2 double or n x 1 cell array of m x 2 double. MS2 spectrum/spectra.
% fragmodes: string or n x 1 cell array of strings. Fragmentation method of spectrum/spectra.
% charges: 1 x 1 double or n x 1 double. Charge state of precursor ion(s).
% monomasses: 1 x 1 double or n x 1 double. Monoisotopic mass of precursor ion(s).
% ms2tols: 1 x 1 double or n x 1 double. MS2 tolerence to be used in peak matching (value).
% ms2tolunits: string or n x 1 cell array of strings. MS2 tolerence to be used in peak matching (Unit).
% denoisingoptions: structure. All scoring parameters. Three fields are used:
%     minmaxmz: 1 x 2 double. The lower and upper limit of m/z of the peaks
%         to be kept. Peaks with m/z out of this range will be removed from
%         output.
%     fragmode: 1 x n cell array of strings. Fragmentation modes to be
%         expected. Each fragmentation mode has its own "cutoffmed"
%         and "fracmax" values.
%     cutoffmed: 1 x n double. Median of peak intensities will be multiplied by this
%         factor, the product is one of the noise filtering thresholds.
%     fracmax: 1 x n double. Highest peak intensity will be multiplied by this
%         factor, the product is one of the noise filtering thresholds.
%
% Output:
% spectra: m x 2 double or n x 1 cell array of m x 2 double. Spectrum/spectra
%     with noise removed.
%
% Note:
% Input can be for 1 spectrum or multiple spectra. In the case of multiple
%     spectra, each spectrum must be stored in individual cells. The
%     corresponding information, including "fragmodes", "charges",
%     "monomasses", "ms22tols", "ms2tolunits" must be prepared for each
%     spectrum separately.
% "cutoffmed" and "fracmax" provides two thresholds. The lower one will
%     be adopted as the working threshold, all peaks whose intensity are
%     smaller than this threshold will be deleted.
%
% Example:
% N/A
%
% Children function:
% REMOVEPRECURSORION  POLISHSPECTRA
%

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 11/10/2020

cutoffmed = denoisingoptions.cutoffmed;
fracmax = denoisingoptions.fracmax;
fragmode = denoisingoptions.fragmode;
if ischar(fragmode)
    fragmode = {fragmode};
end
if isfield(denoisingoptions,'minmaxmz')
    minmaxmz = denoisingoptions.minmaxmz;
else
    minmaxmz = [-Inf,Inf];
end

% If program is treating multiple spectra
if iscell(spectra)
    for ii = 1:length(spectra)
        spectrum = spectra{ii};
        tempfragmode = datafragmodes{ii};
        if ~isempty(spectrum)
            % REMOVEPRECURSORION is applied on ETD-based spectra only
            % Remove precursor ions
            if ismember(upper(tempfragmode),{'ETD','ETCID','ETHCD'})
                tempprecmh = datamonomasses(ii) + 1.007825032;
                spectrum = removePrecursorIon(spectrum,tempprecmh,datacharges(ii),ms2tols(ii),...
                    ms2tolunits{ii});
            end
            tempcutoffmed = cutoffmed(ismember(upper(fragmode),upper(tempfragmode)));
            tempfracmax = fracmax(ismember(upper(fragmode),upper(tempfragmode)));
            % Remove noise peaks
            spectrum = PolishSpectra(spectrum,tempcutoffmed,tempfracmax);
            % Apply m/z lower and upper limits
            keepthispeak = (spectrum(:,1) >= minmaxmz(1)) & (spectrum(:,1) <= minmaxmz(2));
            spectra{ii} = spectrum(keepthispeak,:);
        end
    end
    
    % If program is treating a single spectrum
elseif isnumeric(spectra)
    if ~isempty(spectra)
        % REMOVEPRECURSORION is applied on ETD-based spectrum only
        % Remove precursor ions
        if ismember(upper(datafragmodes),{'ETD','ETCID','ETHCD'})
            tempprecmh = datamonomasses + 1.007825032;
            spectra = removePrecursorIon(spectra,tempprecmh,datacharges,ms2tols,...
                ms2tolunits);
        end
        if (length(cutoffmed)) > 1 && (length(fracmax) > 1)
            tempcutoffmed = cutoffmed(ismember(upper(fragmode),upper(datafragmodes)));
            tempfracmax = fracmax(ismember(upper(fragmode),upper(datafragmodes)));
        else
            tempcutoffmed = cutoffmed;
            tempfracmax = fracmax;
        end
            
        % Remove noise peaks
        spectra = PolishSpectra(spectra,tempcutoffmed,tempfracmax);
        % Apply m/z lower and upper limits
        keepthispeak = (spectra(:,1) >= minmaxmz(1)) & (spectra(:,1) <= minmaxmz(2));
        spectra = spectra(keepthispeak,:);
    end
end
end