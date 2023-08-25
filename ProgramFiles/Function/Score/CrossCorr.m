function [peakLag,htCenter,htAvg]=CrossCorr(SpectraA,SpectraB,maxlag,tolUnit)
% CROSSCORR: Perform a cross-correlation analysis between SpectraA and SpectraB
% to determine how well they are correlated
%
% Syntax:
%    [peakLag,htCenter,htAvg]=CrossCorr(SpectraA,SpectraB,maxlag,tolUnit)
%
% Input: 'Spectra A' and 'Spectra B' which have to be compared, and
% maxLag, the distance by which they are offset initially.
%
% Output: Peak value of the normalized Crosscorrelation function at the center (i.e
% htCenter at lag <= |1|), ht at center with respect to the overall mean, and the
% offset when this function peaks.
%
% Children Function: None
%
% See also scoreAllSpectra, Pvalue, scoreProb, score1Spectra,
% scoreAllSpectra.

% Author: Sriram Neelamegham
% Date Lastly Updated: 8/11/14 by Gang Liu

% Note: For more infor regarding xcorr, see:
% "http://www.mathworks.com/matlabcentral/answers/19652-xcorr"
% maxlag=50        % by default, we look for lag=+/-50
% Last edited on 12/22/2013

maxA=ceil(max(SpectraA(:,1)));  % find the range of SpectraA
if strcmpi(tolUnit,'Da')
    mode=1;
else
    mode=2;
end

switch mode
    case 1   % for lower resolution when units are Da
        A=zeros(maxA+1,1);  % initialize A and B, arrays with one element for each m/z
        B=A;
        
        for i=1:size(SpectraA,1)
            bottom=floor(SpectraA(i,1));         % even if high resolution we digitize by 1Da
            top=ceil(SpectraA(i,1));
            if (bottom>0)
                A(bottom)=A(bottom)+SpectraA(i,2); % This adds Spectra value to upper and lower bounds of A
                A(top)=A(top)+SpectraA(i,2);       % Note: if top=bottom, this is added twice and this is ok
            end
        end
        
        for i=1:size(SpectraB,1)
            bottom = floor(SpectraB(i,1));
            top    = ceil(SpectraB(i,1));
            if ((bottom>1) && (top <= maxA))
                B(bottom) = B(bottom)+SpectraB(i,2);    % This adds Spectra value to upper and lower bounds of B
                B(top)    = B(top)+SpectraB(i,2);       % Note: if top=bottom, this is added twice and this is ok
            end
        end
        
    case 2    % for higher resolution when units are ppm
        A=zeros(maxA*20+1,1);  % initialize A and B, arrays with one element for each m/z
        B=A;
        
        for i=1:size(SpectraA,1)
            bottom=floor(SpectraA(i,1)*20);
            top=ceil(SpectraA(i,1)*20);
            if (bottom>0)
                A(bottom)=A(bottom)+SpectraA(i,2); % This adds Spectra value to upper and lower bounds of A
                A(top)=A(top)+SpectraA(i,2);       % Note: if top=bottom, this is added twice and this is ok
            end
        end
        
        for i=1:size(SpectraB,1)
            bottom = floor(SpectraB(i,1)*20);
            top    = ceil(SpectraB(i,1)*20);
            if ((bottom>1) && (top <= maxA*20))
                B(bottom) = B(bottom)+SpectraB(i,2);    % This adds Spectra value to upper and lower bounds of B
                B(top)    = B(top)+SpectraB(i,2);       % Note: if top=bottom, this is added twice and this is ok
            end
        end
        if length(A) > 8001
            A = A(8001:end,:);  % 400 Da cut-off
            B = B(8001:end,:);
        end
end
[c,lags] = xcorr(A,B,maxlag,'coeff');   % option 'co-eff' provides normalized plot where peak Xcorr=1
if any(isnan(c))
    peakLag = -50;
    htCenter = 0;
    htAvg    = 0;
else
    peakLag  = lags(c==max(c));          % Lag value at which XCorr function peaks
    peakLag  = max(abs(peakLag));
    htCenter = max(c(maxlag:maxlag+2));       % ht. at center
    htAvg    = htCenter/mean(c);                 % ht. at center with respect to overall mean
end
end