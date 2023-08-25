function [spectrascores,stopsignal] = presearch_match1byn(sgp,sgpmass,theofrag,proteinID,scannums,...
    exptmass,monomass,retime,spectra,charge,...
    tol,tolunit,fragmethod,fragnum,mdiff,...
    quant,proteinname,options)
% PRESEARCH_MATCH1BY1: get the score for 1 spectrum vs 1 candidate for 1
% fragmentation mode
%
% Syntax:
% [spectrascores,stopsignal] = presearch_match1by1(sgp,sgpmass,theofrag,proteinID,...
%     scannum,exptmass,monomass,retime,spectrum,charge,tol,tolunit,fragmethod,...
%     fragnum,options,stopsignal)
%
% Input:
% sgp: string. SGP sequence of candidate (glyco)peptide.
% sgpmass: mass of this SGP (not M+H)
% theofrag: structure. Theoretical fragments of candidate (glyco)peptide.
% proteinID: 1 x n numerical array. Indices of candidate.
%     [Protein #, Peptide #, PTM #1, PTM#2,...]
% scannum: Integer. Spectrum's scan number.
% exptmass: mass of the ion that got fragmented.
% monomass: (proposed) monoisotopic mass of the ion.
% retime: retention time of this scan (MS2).
% spectrum: n x 2 numerical array. Spectrum to be analyzed.
% charge: Integer. Charge carried by parent ion.
% tol: non-zero number. Tolerence value.
% tolunit: string. Tolerence unit, "Da" or "ppm".
% fragmethod: string. Fragmentation method of parent ion.
% fragnum: [NpFrag, NgFrag, NmFrag]
% options: structure with fields:
%     maxlag: maximum lag for cross correlation calculation
%     selectpeak: additional m/z's to search
% stopsignal: logical value. When set to "true", function will return a
% default empty result.
%
% Output:
% spectrascores: structure. Scores from this match. Fields are:
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%     Field name                        Format                                 Description
%     scannum                            n x 1 double                        Scan number of spectra
%     Scan                                    n x 1 double                        Scan number of spectrum
%     Mono                                  n x 1 double                       Monoisotopic (experimental) mass of precursor ion
%     Theo                                    n x 1 double                       Theoretical mass of candidate
%     Expt                                     n x 1 double                       Actual mass of precursor ion that was fragmented
%     Charge                               n x 1 double                        Charge state of precursor ion
%     SGP                                     n x 1 cell array of               SGP 1.0 sequence of candidate
%                                                      strings
%     Fragmode                         n x 1 cell array of               Fragmentation method
%                                                      strings
%     PeakLag                             n x 1 double                        Offset between theoretical and experimental spectrum where
%                                                                                                     maximum correlation occur
%     HtCenter                           n x 1 double                        Value of maximum correlationm
%     HtAvg                                 n x 1 double                       "HtCenter" divided by mean of all cross correlation values
%     PercentIonMatch           n x 1 double                       Percentage of theoretical fragments that are matched
%     Pvalue                                n x 1 double                       P-value of candidate
%     DecoyRatio                       n x 1 double                       T-test value of candidate
%     Top10                                n x 1 double                        Number of 10 highest peaks in spectrum that are matched
%     SelectPeak                        n x 1 cell array of              User defined m/z values that are matched
%                                                      strings
%     NpFrag                              n x 1 double                        Number of peptide fragmentation performed
%     NgFrag                              n x 1 double                        Number of glycan fragmentation performedm
%     NmFrag                             n x 1 double                        Number of non-glycan modification fragmentation performed
%     Enscore                             n x 1 double                        Ensemble score of the candidate
%     DecoyEnscore                 n x 1 double                        Ensemble score of the decoy of the candidate
%     Retime                              n x 1 double                         Retention time of scan
%     Quant                                n x 1 double                        Result of label free quantification
%     Glycov                               n x 1 double                        Percentage of glycosidic bond that are fragmented
%     Y0Y1Y2                             n x 1 cell array of                         m
%                                                     strings
%     ProteinID                         n x 1 double m
%     Fracpepfragmatch        n x 1 double m
%     Fracglyfragmatch          n x 1 double m
%     MDiff                                n x 1 double m
%     Protein                             n x 1 cell array of s
%                                                    strings
%
%
%
%
%
% ---------------------------------------------------------------------------------------------------------------------------------------------------
% stopsignal: logical value. Originally designed for "trigger" experiements.
% If during analysis this candidate seems to unlikely, MATCH1BY1 will
% return "true", therefore the rest of the experiement in this trigger
% group may be avoided.
%
% Note:
% It is recommended that spectrum be centroided and contains no zero
%     intensities.
% If the matching result is too poor (Top 10 < 2 and less than 1% of all
% theoretical ions were matched), program will terminate and return the
% default empty result.
%
% Example:
% N/A
%
% Children function:
% CALCITHSCORE  CALCULATETOP10MATCH  CALCPEPCLEAVAGECOV
% CALCY0Y1Y2  CALCGLYCLEAVAGECOV  IONDECOY  COMPCIDENSCORE
% COMPHCDENSCORE  COMPETDENSCORE  COMPETCIDENSCORE
% COMPETHCDENSCORE

% Notes: Line 145 has criteria that may significantly speed up calculations
% if adjusted

%% INITIALIZE OUTPUT (default output that is returned if the match is not good)
stopsignal = true;
cisoptions.maxlag = options.maxlag;
cisoptions.selectpeak = options.selectpeak;
%% CALC START
goodmatch = false(size(theofrag));
for ii = 1:length(theofrag)
    theofragmz = [theofrag{ii}.mz];  % theoretical fragmentation spectrum
    spectrum = spectra{ii};
    thistheores = calcithscore(spectrum,theofragmz,charge,tol(ii),tolunit{ii},2,cisoptions);
    thisTop10 = calculatetop10match(spectrum,thistheores.peakmatchindex,2,400);
    if thisTop10 >= 3
        goodmatch(ii) = true;
    end
end
if ~any(~goodmatch)
    nDecoy = 25;
    stopsignal = false;
    for ii = length(theofrag):-1:1
        spectrum = spectra{ii};
        temptheo
        spectrascores(ii).Scan = scannums(ii);
        spectrascores(ii).Mono = monomass(ii);  % experiment monoisotopic mass
        spectrascores(ii).Theo = sgpmass(ii);  % candidate mass
        spectrascores(ii).Expt = exptmass(ii);  % experiment fragmented mass
        spectrascores(ii).Charge = charge(ii);
        spectrascores(ii).SGP = sgp;
        spectrascores(ii).Fragmode = fragmethod;
        spectrascores(ii).PeakLag  = -50;
        spectrascores(ii).HtCenter = 0;
        spectrascores(ii).HtAvg    = 0;
        spectrascores(ii).PercentIonMatch = 0;
        spectrascores(ii).Pvalue   = 0.99;
        spectrascores(ii).DecoyRatio = 0.99;
        spectrascores(ii).Top10 = 0;
        spectrascores(ii).SelectPeak = '';
        spectrascores(ii).NpFrag = fragnum(1);
        spectrascores(ii).NgFrag = fragnum(2);
        spectrascores(ii).NmFrag = fragnum(3);
        spectrascores(ii).Enscore = 0;
        spectrascores(ii).DecoyEnscore = 0;
        spectrascores(ii).Retime = retime(ii);
        spectrascores(ii).Quant = quant(ii);
        spectrascores(ii).Glycov = 0;
        spectrascores(ii).Pepcov = 0;
        spectrascores(ii).Y0Y1Y2 = num2str([0,0,0]);
        spectrascores(ii).ProteinID = proteinID;
        spectrascores(ii).Fracpepfragmatch = 0;
        spectrascores(ii).Fracglyfragmatch = 0;
        spectrascores(ii).MDiff = mdiff;
        spectrascores(ii).Protein = proteinname;
        
        theofragmz = [theofrag{ii}.mz];  % theoretical fragmentation spectrum
        theores = calcithscore(spectrum,theofragmz,charge,tol(ii),tolunit{ii},1,cisoptions);
        [PGM{1},PGM{2},PGM{3}] = breakGlyPep(theofrag(1).original);  % these are pepMat,glyMat,modMat
        spectrascores(ii).SelectPeak = num2str(cisoptions.selectpeak(theores.selectpeakismatched));
        fragAAind = cell(size(theofrag{ii}));
        fragglyind = cell(size(theofrag{ii}));
        fragmodind = cell(size(theofrag{ii}));
        for jj = 1:length(fragAAind)
            fragAAind{jj} = theofrag(jj).unitindex{1};
            fragglyind{jj} = theofrag(jj).unitindex{2};
            fragmodind{jj} = theofrag(jj).unitindex{3};
        end
        theofragtype = {theofrag{ii}.type};
        %% PART - PEPTIC BOND CLEAVAGE CALCULATION (w/wo GLYCAN)
        peplen = length(PGM{1}.pep);
        if ~isempty(PGM{2})
            glypos = [PGM{2}.pos];
        else
            glypos = 0;
        end
        [pfismatched,gfismatched,pepcov] = calcpepcleavagecov(theores,theofrag,fragAAind,peplen,glypos);
        spectrascores(ii).Fracpepfragmatch = sum(pfismatched) / (peplen - 1);
        spectrascores(ii).Fracglyfragmatch = sum(gfismatched) / (peplen - 1);
        spectrascores(ii).Pepcov = sum(logical(pfismatched|gfismatched))/(peplen - 1);
        %% PART - Y0Y1Y2
        Y0Y1Y2 = calcY0Y1Y2(spectrum,theores,fragglyind,fragAAind,peplen);  % Y0Y1Y2 CALCULATION
        %% GLYCAN
        glycov = 1;
        if ~isempty(PGM{2})
            glycov = calcglycleavagecov(PGM,fragglyind,theores);
            spectrascores(ii).Glycov = glycov;
        end
        %% PART - OXO
        isoxofrag = ~cellfun(@isempty,strfind(theofragtype,'oxo')) | ismember(theofragtype,{'-b{n}','-b{h}','-b{s}','-b{n{h}}'});
        matchedoxo = reshape(~cellfun(@isempty,theores.ionmatchindex),[],1) & reshape(isoxofrag,[],1);
        percentoxo = sum(matchedoxo)/sum(isoxofrag);
        if isnan(percentoxo)
            percentoxo = 0;
        end
        %% PART - DECOY
        decoyfound = 0;
        decoysearch = 0;
        decoyArray = zeros(1,nDecoy);
        for jj = 1:nDecoy
            decoyfragmz = Iondecoy(PGM,'RANDOM',theofragmz,fragAAind,fragglyind,fragmodind,'glypep');
            decoyres = calcithscore(spectrum,decoyfragmz,charge,tol(ii),tolunit{ii},3,cisoptions);
            ionismatched = decoyres.ionmatchindex;
            found = sum(ionismatched);
            search = length(ionismatched);
            decoyfound = decoyfound + found;
            decoysearch = decoysearch + search;
            decoyArray(jj) = found/search * 100;
        end
        decoyres = calcithscore(spectrum,decoyfragmz,charge,tol(ii),tolunit{ii},1,cisoptions);
        %% calculate score
        % P-value
        pnorm = decoyfound/decoysearch;
        if theores.found == 0
            Pvalue = 1;
        else
            Pvalue = 1 - poisscdf(theores.found,theores.search*pnorm);
        end
        if decoyres.found == 0
            decoyPvalue = 1;
        else
            decoyPvalue = 1 - poisscdf(decoyres.found,decoyres.search*pnorm);
        end
        decoyRatio = ttest2(decoyres.percentIonMatch,decoyArray);
        if isnan(decoyRatio)
            decoyRatio = 0.99;
        end
        decoyglycov = 1;
        if ~isempty(PGM{2})
            decoyglycov = calcglycleavagecov(PGM,fragglyind,decoyres);
        end
        [~,~,decoypepcov] = calcpepcleavagecov(decoyres,theofrag,fragAAind,peplen,glypos);
        decoyY0Y1Y2 = calcY0Y1Y2(spectrum,decoyres,fragglyind,fragAAind,peplen);
        decoymatchedoxo = reshape(~cellfun(@isempty,decoyres.ionmatchindex),[],1) & reshape(isoxofrag,[],1);
        decoypercentoxo = sum(decoymatchedoxo)/sum(isoxofrag);
        if isnan(decoypercentoxo)
            decoypercentoxo = 0;
        end
        
        if strcmpi(fragmethod,'CID')
            enscore = compCIDenscore(glycov,theores.peakLag,Pvalue,theores.htAvg,theores.Top10);
            decoy_enscore = compCIDenscore(decoyglycov,decoyres.peakLag,decoyPvalue,decoyres.htAvg,decoyres.Top10);
        elseif strcmpi(fragmethod,'HCD') || strcmpi(fragmethod,'HCDfragilemonosac')
            enscore = compHCDenscore(pepcov,theores.peakLag,Pvalue,theores.Top10,Y0Y1Y2,percentoxo);
            decoy_enscore = compHCDenscore(decoypepcov,decoyres.peakLag,decoyPvalue,decoyres.Top10,decoyY0Y1Y2,decoypercentoxo);
        elseif strcmpi(fragmethod,'ETD')
            enscore = compETDenscore(pepcov,theores.peakLag,theores.htAvg,Pvalue);
            decoy_enscore = compETDenscore(decoypepcov,decoyres.peakLag,decoyres.htAvg,decoyPvalue);
        elseif strcmpi(fragmethod,'ETCID')
            % enscore = compETCIDenscore(theores.peakLag,theores.htAvg,Pvalue,theores.Top10,theores.percentIonMatch);
            % decoy_enscore = compETCIDenscore(decoyres.peakLag,decoyres.htAvg,decoyPvalue,decoyres.Top10,decoyres.percentIonMatch);
            enscore = compETHCDenscore(glycov,pepcov,theores.peakLag,Pvalue,theores.htCenter,theores.Top10,percentoxo);
            decoy_enscore = compETHCDenscore(decoyglycov,decoypepcov,decoyres.peakLag,decoyPvalue,decoyres.htCenter,decoyres.Top10,decoypercentoxo);
        elseif strcmpi(fragmethod,'ETHCD')
            enscore = compETHCDenscore(glycov,pepcov,theores.peakLag,Pvalue,theores.htCenter,theores.Top10,percentoxo);
            decoy_enscore = compETHCDenscore(decoyglycov,decoypepcov,decoyres.peakLag,decoyPvalue,decoyres.htCenter,decoyres.Top10,decoypercentoxo);
        end
        spectrascores(ii).PeakLag = theores.peakLag;
        spectrascores(ii).HtCenter = theores.htCenter;
        spectrascores(ii).HtAvg = theores.htAvg;
        spectrascores(ii).Pvalue = Pvalue;
        spectrascores(ii).PercentIonMatch = theores.percentIonMatch;
        spectrascores(ii).DecoyRatio = decoyRatio;
        spectrascores(ii).Top10 = theores.Top10;
        spectrascores(ii).Enscore = enscore;
        spectrascores(ii).DecoyEnscore = decoy_enscore;
        spectrascores(ii).Y0Y1Y2 = num2str(Y0Y1Y2);
    end
end
end

function enscore=compCIDenscore(glycov,peakLag,Pvalue,htAvg,Top10)
%
% calculate Ensemble Score for CID mode
%
if(abs(peakLag)<=1)
    if htAvg > 10
        norm_xcorr=1;
    else
        norm_xcorr=min(htAvg/10,1);
    end
else
    norm_xcorr=0;
end

pvaluethreshold=1e-5;
if((Pvalue-0.02)>0)
    pscore=0;
elseif(Pvalue<pvaluethreshold)
    pscore=1;
else
    kslope = 1/(log(pvaluethreshold)-log(0.02));
    pscore = 1+kslope*(log(Pvalue)-log(pvaluethreshold));
end

if glycov > .8
    perglycov = 1;
elseif glycov < .2
    perglycov = 0;
else
    perglycov = (glycov-.2)/.6;
end

if Top10 >= 7
    Top10score = 1;
else
    Top10score = Top10/7;
end

top10scoreweight = 0.25; pscoreweight=0.25;
glycovweight = 0.25; xcorrweight = 0.25;
enscore = Top10score*top10scoreweight + pscore*pscoreweight + ...
    perglycov*glycovweight + norm_xcorr*xcorrweight;
if(isnan(enscore))
    enscore =0;
end

end

function enscore=compHCDenscore(pepcov,peakLag,Pvalue,Top10,Y012,Oxocov)
%
% calculate Ensemble Score for HCD mode
%
pvaluethreshold=1e-5;
if((Pvalue-0.02)>0)
    pscore=0;
elseif(Pvalue<pvaluethreshold)
    pscore=1;
else
    kslope = 1/(log(pvaluethreshold)-log(0.02));
    pscore = 1+kslope*(log(Pvalue)-log(pvaluethreshold));
end

if(sum(Y012) > .15)
    Y012match = 1;
else
    Y012match = sum(Y012)/.15;
end

top10scoreweight = 1/6; pscoreweight=1/6;
Y012weight = 1/3; pepcovweight = 1/6;
Oxoweight = 1/6;

if Top10 >= 7
    Top10score = 1;
else
    Top10score = Top10/7;
end

enscore = Top10score*top10scoreweight + pscore*pscoreweight + ...
    Y012match*Y012weight + pepcov*pepcovweight + Oxocov*Oxoweight;
if(isnan(enscore))
    enscore =0;
end

end

function enscore=compETDenscore(pepcov,peakLag,htCenter,Pvalue)
%
% calculate Ensemble Score for ETD mode
%

htcenternorm = 0.50;
if(abs(peakLag)<=1)
    if htCenter > htcenternorm
        norm_xcorr = 1;
    else
        norm_xcorr = (htCenter/htcenternorm);
    end
else
    norm_xcorr = 0;
end

pvaluethreshold=1e-5;
pvalueupper = 0.1;
if((Pvalue - pvalueupper)>0)
    pscore = 0;
elseif(Pvalue<pvaluethreshold)
    pscore = 1;
else
    kslope = 1/(log(pvaluethreshold)-log(pvalueupper));
    pscore = 1+kslope*(log(Pvalue)-log(pvaluethreshold));
end

pepcovscoreweight  = 1/3; pscoreweight = 1/3;
xcorrweight  = 1/3;

enscore = pepcov*pepcovscoreweight + pscore*pscoreweight + ...
    norm_xcorr*xcorrweight;

if(isnan(enscore))
    enscore =0;
end
if(enscore<0)
    enscore = 0;
end

end

function enscore=compETCIDenscore(peakLag,htCenter,Pvalue,Top10,percentIonMatch)
%
% calculate Ensemble Score for ETD mode
%

htcenternorm = 0.50;
if(abs(peakLag)<=1)
    if htCenter > htcenternorm
        norm_xcorr = 1;
    else
        norm_xcorr = (htCenter/htcenternorm);
    end
else
    norm_xcorr = 0;
end

pvaluethreshold=1e-5;
pvalueupper = 0.1;
if((Pvalue - pvalueupper)>0)
    pscore = 0;
elseif(Pvalue<pvaluethreshold)
    pscore = 1;
else
    kslope = 1/(log(pvaluethreshold)-log(pvalueupper));
    pscore = 1+kslope*(log(Pvalue)-log(pvaluethreshold));
end

if(percentIonMatch>50)
    perIonmatch=1;
else
    perIonmatch=percentIonMatch/50;
end

if Top10 >= 5
    Top10score = 1;
else
    Top10score = Top10/5;
end

top10scoreweight  = 0.10; pscoreweight = 0.70;
perIonmatchweight = 0;    xcorrweight  = 0.20;

enscore = Top10score*top10scoreweight + pscore*pscoreweight + ...
    perIonmatch*perIonmatchweight + norm_xcorr*xcorrweight;

if(isnan(enscore))
    enscore =0;
end

if(enscore<0)
    enscore = 0;
end
end

function enscore=compETHCDenscore(glycov,pepcov,peakLag,Pvalue,htCenter,Top10,Oxocov)
% peakLag,htCenter,Pvalue,Top10,percentIonMatch,nummatchedpeaks)
%
% calculate Ensemble Score for ETD mode
%

htcenternorm = 0.50;
if(abs(peakLag)<=1)
    if htCenter > htcenternorm
        norm_xcorr = 1;
    else
        norm_xcorr = (htCenter/htcenternorm);
    end
else
    norm_xcorr = 0;
end

pvaluethreshold=1e-5;
pvalueupper = 0.1;
if((Pvalue - pvalueupper)>0)
    pscore = 0;
elseif(Pvalue<pvaluethreshold)
    pscore = 1;
else
    kslope = 1/(log(pvaluethreshold)-log(pvalueupper));
    pscore = 1 + kslope*(log(Pvalue)-log(pvaluethreshold));
end

glycovweight = 1/8;
pepcovweight = 1/4;
pscoreweight = 1/8;
top10scoreweight = 1/4;
xcorrweight  = 1/8;
Oxocovweight = 1/8;

if Top10 >= 5
    Top10score = 1;
else
    Top10score = Top10/5;
end

enscore = glycov*glycovweight + pepcov*pepcovweight + ...
    pscore*pscoreweight + norm_xcorr*xcorrweight + ...
    Top10score*top10scoreweight + Oxocov*Oxocovweight;

if(isnan(enscore))
    enscore =0;
end
if(enscore<0)
    enscore = 0;
end
end