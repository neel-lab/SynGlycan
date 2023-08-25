function specsimp = ms2centroid(specin,tol,tolUnit,premass)
%MS2CENTROID: Combines peaks that are within the mass tolerance in MS2 and
% MS3 spectra
% 
% Syntax:
%   specsimp = ms2centroid(specin)
% 
% Input:
%   specin: n x 2 numerical array, input spectrum
% 
% Output:
%   specsimp: n x 2 numerical array, centroided spectrum
% 
% Note: Mass tolerance must be changed directly in the code
%
%See also: QUANTITATIVEANALYSIS
%

specin_og = specin;
try
    %tol = 1; % mass tolerance
    %tolUnit = 'Da'; % mass tolerance units
    specold = specin; % input spectra (m/z and intensity values)
    specnew = [0,0]; % index value
        if strcmpi(tolUnit,'ppm')  % for ppm
             tolC = tol/1e6*premass;  % mass tolerance in ppm
        elseif (strcmpi(tolUnit,'Da')) % for Da
             tolC = tol;  % mass tolerance in Da   
        else    
             error('MATLAB:GLYCOPAT:ERRORUNIT','INCORRECT UNIT FOR MASS TOL');
        end 
        while length(specold) ~= length(specnew)
            specold = specin;
            for ii = 2:length(specin)
                jj = ii-1;
                if (specin(ii,1) - specin(jj,1)) < tolC
                    e = max(specin(jj,2),specin(ii,2));
                    f = find (specin(:,2) == e);
                    if f == ii
                        specin(ii,2) = specin(ii,2) + specin(jj,2);
                        specin(jj,:) = [0,0];
                    else
                        specin(jj,2) = specin(ii,2) + specin(jj,2);
                        specin(ii,:) = [0,0];
                    end
                end
            end
            specin = specin(any(specin,2),:);
            specnew = specin;
         end    
         specsimp = specin;
catch
    specsimp = specin_og;
end



















