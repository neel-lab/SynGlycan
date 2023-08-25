function [spectrascores,stopsignal] = match1by1(sgp,sgpmass,theofrag,proteinID,scannum,...
    exptmass,monomass,retime,spectrum,charge,...
    tol,tolunit,fragmethod,fragnum,mdiff,...
    quant,proteinname,options,stopsignal)
% MATCH1BY1: get the score for 1 spectrum vs 1 candidate for 1
% fragmentation mode
%
% Syntax:
% [spectrascores,stopsignal] = match1by1(sgp,sgpmass,theofrag,proteinID,...
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
spectrascores.Scan = scannum;
spectrascores.Mono = monomass;  % experiment monoisotopic mass
spectrascores.Theo = sgpmass;  % candidate mass
spectrascores.Expt = exptmass;  % experiment fragmented mass
spectrascores.Charge = charge;
spectrascores.SGP = sgp;
spectrascores.Fragmode = fragmethod;
spectrascores.PeakLag  = -50;
spectrascores.HtCenter = 0;
spectrascores.HtAvg    = 0;
spectrascores.PercentIonMatch = 0;
spectrascores.Pvalue   = 0.99;
spectrascores.DecoyRatio = 0.99;
spectrascores.Top10 = 0;
spectrascores.SelectPeak = '';
spectrascores.NpFrag = fragnum(1);
spectrascores.NgFrag = fragnum(2);
spectrascores.NmFrag = fragnum(3);
spectrascores.Enscore = 0;
spectrascores.DecoyEnscore = 0;
spectrascores.Retime = retime;
spectrascores.Quant = quant;
spectrascores.Glycov = 0;
spectrascores.Pepcov = 0;
spectrascores.Y0Y1Y2 = num2str([0,0,0]);
spectrascores.ProteinID = proteinID;
spectrascores.Fracpepfragmatch = 0;
spectrascores.Fracglyfragmatch = 0;
spectrascores.MDiff = mdiff;
spectrascores.Protein = proteinname;

nDecoy = GlycoPATConstant.ndecoy;

%% CALC START
cisoptions.maxlag = options.maxlag;
cisoptions.selectpeak = options.selectpeak;
proceed_override = options.proceed_override;
proceed = true;
stopsignal = false;
if size(spectrum,1) >= 5
    theofragmz = [theofrag.mz];  % theoretical fragmentation spectrum
    theores = calcithscore(spectrum,theofragmz,charge,tol,tolunit,1,cisoptions);
    [PGM{1},PGM{2},PGM{3}] = breakGlyPep(theofrag(1).original);  % these are pepMat,glyMat,modMat
    peplen = length(PGM{1}.pep);
    fragAAind = cell(size(theofrag));
    fragglyind = cell(size(theofrag));
    fragmodind = cell(size(theofrag));
    for ii = 1:length(fragAAind)
        fragAAind{ii} = theofrag(ii).unitindex{1};
        fragglyind{ii} = theofrag(ii).unitindex{2};
        fragmodind{ii} = theofrag(ii).unitindex{3};
    end
    %% PART - Y0Y1Y2
    Y0Y1Y2 = calcY0Y1Y2(spectrum,theores,fragglyind,fragAAind,peplen);  % Y0Y1Y2 CALCULATION
    if ~proceed_override && ~any(Y0Y1Y2) && strcmpi(fragmethod,'HCD')
        proceed = false;
        stopsignal = true;
    end
    if proceed
        spectrascores.SelectPeak = num2str(cisoptions.selectpeak(theores.selectpeakismatched));
        theofragtype = {theofrag.type};
        %% PART - PEPTIC BOND CLEAVAGE CALCULATION (w/wo GLYCAN)
        if ~isempty(PGM{2})
            glypos = [PGM{2}.pos];
        else
            glypos = 0;
        end
        [pfismatched,gfismatched,pepcov] = calcpepcleavagecov(theores,theofrag,fragAAind,peplen,glypos);
        spectrascores.Fracpepfragmatch = sum(pfismatched) / (peplen - 1);
        spectrascores.Fracglyfragmatch = sum(gfismatched) / (peplen - 1);
        spectrascores.Pepcov = sum(logical(pfismatched|gfismatched))/(peplen - 1);
        %% GLYCAN
        glycov = 1;
        if ~isempty(PGM{2})
            glycov = calcglycleavagecov(PGM,fragglyind,theores);
            spectrascores.Glycov = glycov;
        end
        %% PART - OXO
        isoxofrag = ~cellfun(@isempty,strfind(theofragtype,'oxo')) | ismember(theofragtype,{'-b{n}','-b{h}','-b{s}','-b{n{h}}'});
        matchedoxo = reshape(~cellfun(@isempty,theores.ionmatchindex),[],1) & reshape(isoxofrag,[],1);
        percentoxo = sum(matchedoxo)/sum(isoxofrag);
        if isnan(percentoxo)
            percentoxo = 0;
        end
        
        %% PART - DECOY
        decoyneeded = true;
        switch upper(fragmethod)
            case 'HCD'
                if GlycoPATConstant.Scoreweight_HCD_Pscore == 0
                    decoyneeded = false;
                end
            case 'CID'
                if GlycoPATConstant.Scoreweight_CID_Pscore == 0
                    decoyneeded = false;
                end
            case 'ETD'
                if GlycoPATConstant.Scoreweight_ETD_Pscore == 0
                    decoyneeded = false;
                end
            case 'ETCID'
                if GlycoPATConstant.Scoreweight_ETciD_Pscore == 0
                    decoyneeded = false;
                end
            case 'ETHCD'
                if GlycoPATConstant.Scoreweight_EThcD_Pscore == 0
                    decoyneeded = false;
                end
            otherwise
                error('MATCH1BY1 - Unsupported Fragmentation Method');
        end
        Pvalue = 0;
        decoyPvalue = 1;
        decoyfragmz = Iondecoy(PGM,'RANDOM',theofragmz,fragAAind,fragglyind,fragmodind,'glypep');
        decoyres = calcithscore(spectrum,decoyfragmz,charge,tol,tolunit,1,cisoptions);
        if decoyneeded
            decoyfound = 0;
            decoysearch = 0;
            decoyArray = zeros(1,nDecoy);
            for ii = 1:nDecoy
                decoyfragmz = Iondecoy(PGM,'RANDOM',theofragmz,fragAAind,fragglyind,fragmodind,'glypep');
                decoyres_temp = calcithscore(spectrum,decoyfragmz,charge,tol,tolunit,3,cisoptions);
                ionismatched = decoyres_temp.ionmatchindex;
                found = sum(ionismatched);
                search = length(ionismatched);
                decoyfound = decoyfound + found;
                decoysearch = decoysearch + search;
                decoyArray(ii) = found/search * 100;
            end
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
        end
        if strcmpi(fragmethod,'CID')
            decoyglycov = calcglycleavagecov(PGM,fragglyind,decoyres);
            enscore = compCIDenscore(glycov,theores.peakLag,Pvalue,theores.htAvg,theores.Top10);
            decoy_enscore = compCIDenscore(decoyglycov,decoyres.peakLag,decoyPvalue,decoyres.htAvg,decoyres.Top10);
        elseif strcmpi(fragmethod,'HCD') || strcmpi(fragmethod,'HCDfragilemonosac')
            [~,~,decoypepcov] = calcpepcleavagecov(decoyres,theofrag,fragAAind,peplen,glypos);
            enscore = compHCDenscore(pepcov,theores.peakLag,Pvalue,theores.Top10,Y0Y1Y2,percentoxo);
            decoyY0Y1Y2 = calcY0Y1Y2(spectrum,decoyres,fragglyind,fragAAind,peplen);
            decoymatchedoxo = reshape(~cellfun(@isempty,decoyres.ionmatchindex),[],1) & reshape(isoxofrag,[],1);
            decoypercentoxo = sum(decoymatchedoxo)/sum(isoxofrag);
            if isnan(decoypercentoxo)
                decoypercentoxo = 0;
            end
            decoy_enscore = compHCDenscore(decoypepcov,decoyres.peakLag,decoyPvalue,decoyres.Top10,decoyY0Y1Y2,decoypercentoxo);
        elseif strcmpi(fragmethod,'ETD')
            enscore = compETDenscore(pepcov,theores.peakLag,theores.htAvg,Pvalue);
            [~,~,decoypepcov] = calcpepcleavagecov(decoyres,theofrag,fragAAind,peplen,glypos);
            decoy_enscore = compETDenscore(decoypepcov,decoyres.peakLag,decoyres.htAvg,decoyPvalue);
        elseif strcmpi(fragmethod,'ETCID')
            decoyglycov = calcglycleavagecov(PGM,fragglyind,decoyres);
            [~,~,decoypepcov] = calcpepcleavagecov(decoyres,theofrag,fragAAind,peplen,glypos);
            enscore = compETHCDenscore(glycov,pepcov,theores.peakLag,Pvalue,theores.htCenter,theores.Top10,percentoxo);
            decoymatchedoxo = reshape(~cellfun(@isempty,decoyres.ionmatchindex),[],1) & reshape(isoxofrag,[],1);
            decoypercentoxo = sum(decoymatchedoxo)/sum(isoxofrag);
            if isnan(decoypercentoxo)
                decoypercentoxo = 0;
            end
            decoy_enscore = compETHCDenscore(decoyglycov,decoypepcov,decoyres.peakLag,decoyPvalue,decoyres.htCenter,decoyres.Top10,decoypercentoxo);
        elseif strcmpi(fragmethod,'ETHCD')
            decoyglycov = calcglycleavagecov(PGM,fragglyind,decoyres);
            [~,~,decoypepcov] = calcpepcleavagecov(decoyres,theofrag,fragAAind,peplen,glypos);
            enscore = compETHCDenscore(glycov,pepcov,theores.peakLag,Pvalue,theores.htCenter,theores.Top10,percentoxo);
            decoymatchedoxo = reshape(~cellfun(@isempty,decoyres.ionmatchindex),[],1) & reshape(isoxofrag,[],1);
            decoypercentoxo = sum(decoymatchedoxo)/sum(isoxofrag);
            if isnan(decoypercentoxo)
                decoypercentoxo = 0;
            end
            decoy_enscore = compETHCDenscore(decoyglycov,decoypepcov,decoyres.peakLag,decoyPvalue,decoyres.htCenter,decoyres.Top10,decoypercentoxo);
        end
        decoyRatio = decoy_enscore/enscore;
        spectrascores.PeakLag = theores.peakLag;
        spectrascores.HtCenter = theores.htCenter;
        spectrascores.HtAvg = theores.htAvg;
        spectrascores.Pvalue = Pvalue;
        spectrascores.PercentIonMatch = theores.percentIonMatch;
        spectrascores.DecoyRatio = decoyRatio;
        spectrascores.Top10 = theores.Top10;
        spectrascores.Enscore = enscore;
        spectrascores.DecoyEnscore = decoy_enscore;
        spectrascores.Y0Y1Y2 = num2str(Y0Y1Y2);
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

% top10scoreweight = 0.25;
% pscoreweight=0.25;
% glycovweight = 0.25;
% xcorrweight = 0.25;
Oxoweight = GlycoPATConstant.Scoreweight_CID_Oxo;
Y012weight = GlycoPATConstant.Scoreweight_CID_Y012;
Glycovweight = GlycoPATConstant.Scoreweight_CID_Glycov;
Pepcovweight = GlycoPATConstant.Scoreweight_CID_Pepcov;
PerIonMatchweight = GlycoPATConstant.Scoreweight_CID_PerIonMatch;
Pscoreweight = GlycoPATConstant.Scoreweight_CID_Pscore;
Top10scoreweight = GlycoPATConstant.Scoreweight_CID_Top10;
XCorrweight = GlycoPATConstant.Scoreweight_CID_XCorr;
sumweight = Oxoweight + Y012weight + Glycovweight + Pepcovweight + PerIonMatchweight + ...
    Pscoreweight + Top10scoreweight + XCorrweight;
Oxoweight = Oxoweight/sumweight;
Y012weight = Y012weight/sumweight;
Glycovweight = Glycovweight/sumweight;
Pepcovweight = Pepcovweight/sumweight;
PerIonMatchweight = PerIonMatchweight/sumweight;
Pscoreweight = Pscoreweight/sumweight;
Top10scoreweight = Top10scoreweight/sumweight;
XCorrweight = XCorrweight/sumweight;

enscore = Top10score*Top10scoreweight + pscore*Pscoreweight + ...
    perglycov*Glycovweight + norm_xcorr*XCorrweight;
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

% top10scoreweight = 1/6; pscoreweight= 0;
% Y012weight = 1/3; pepcovweight = 1/3;
% Oxoweight = 1/6;
Oxoweight = GlycoPATConstant.Scoreweight_HCD_Oxo;
Y012weight = GlycoPATConstant.Scoreweight_HCD_Y012;
Glycovweight = GlycoPATConstant.Scoreweight_HCD_Glycov;
Pepcovweight = GlycoPATConstant.Scoreweight_HCD_Pepcov;
PerIonMatchweight = GlycoPATConstant.Scoreweight_HCD_PerIonMatch;
Pscoreweight = GlycoPATConstant.Scoreweight_HCD_Pscore;
Top10scoreweight = GlycoPATConstant.Scoreweight_HCD_Top10;
XCorrweight = GlycoPATConstant.Scoreweight_HCD_XCorr;
sumweight = Oxoweight + Y012weight + Glycovweight + Pepcovweight + PerIonMatchweight + ...
    Pscoreweight + Top10scoreweight + XCorrweight;
Oxoweight = Oxoweight/sumweight;
Y012weight = Y012weight/sumweight;
Glycovweight = Glycovweight/sumweight;
Pepcovweight = Pepcovweight/sumweight;
PerIonMatchweight = PerIonMatchweight/sumweight;
Pscoreweight = Pscoreweight/sumweight;
Top10scoreweight = Top10scoreweight/sumweight;
XCorrweight = XCorrweight/sumweight;


if Top10 >= 7
    Top10score = 1;
else
    Top10score = Top10/7;
end

enscore = Top10score*Top10scoreweight + pscore*Pscoreweight + ...
    Y012match*Y012weight + pepcov*Pepcovweight + Oxocov*Oxoweight;
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

% pepcovscoreweight  = 1/3; pscoreweight = 1/3;
% xcorrweight  = 1/3;
Oxoweight = GlycoPATConstant.Scoreweight_ETD_Oxo;
Y012weight = GlycoPATConstant.Scoreweight_ETD_Y012;
Glycovweight = GlycoPATConstant.Scoreweight_ETD_Glycov;
Pepcovweight = GlycoPATConstant.Scoreweight_ETD_Pepcov;
PerIonMatchweight = GlycoPATConstant.Scoreweight_ETD_PerIonMatch;
Pscoreweight = GlycoPATConstant.Scoreweight_ETD_Pscore;
Top10scoreweight = GlycoPATConstant.Scoreweight_ETD_Top10;
XCorrweight = GlycoPATConstant.Scoreweight_ETD_XCorr;
sumweight = Oxoweight + Y012weight + Glycovweight + Pepcovweight + PerIonMatchweight + ...
    Pscoreweight + Top10scoreweight + XCorrweight;
Oxoweight = Oxoweight/sumweight;
Y012weight = Y012weight/sumweight;
Glycovweight = Glycovweight/sumweight;
Pepcovweight = Pepcovweight/sumweight;
PerIonMatchweight = PerIonMatchweight/sumweight;
Pscoreweight = Pscoreweight/sumweight;
Top10scoreweight = Top10scoreweight/sumweight;
XCorrweight = XCorrweight/sumweight;


enscore = pepcov*Pepcovweight + pscore*Pscoreweight + ...
    norm_xcorr*XCorrweight;

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

% top10scoreweight  = 0.10; pscoreweight = 0.70;
% perIonmatchweight = 0;    xcorrweight  = 0.20;
Oxoweight = GlycoPATConstant.Scoreweight_ETciD_Oxo;
Y012weight = GlycoPATConstant.Scoreweight_ETciD_Y012;
Glycovweight = GlycoPATConstant.Scoreweight_ETciD_Glycov;
Pepcovweight = GlycoPATConstant.Scoreweight_ETciD_Pepcov;
PerIonMatchweight = GlycoPATConstant.Scoreweight_ETciD_PerIonMatch;
Pscoreweight = GlycoPATConstant.Scoreweight_ETciD_Pscore;
Top10scoreweight = GlycoPATConstant.Scoreweight_ETciD_Top10;
XCorrweight = GlycoPATConstant.Scoreweight_ETciD_XCorr;
sumweight = Oxoweight + Y012weight + Glycovweight + Pepcovweight + PerIonMatchweight + ...
    Pscoreweight + Top10scoreweight + XCorrweight;
Oxoweight = Oxoweight/sumweight;
Y012weight = Y012weight/sumweight;
Glycovweight = Glycovweight/sumweight;
Pepcovweight = Pepcovweight/sumweight;
PerIonMatchweight = PerIonMatchweight/sumweight;
Pscoreweight = Pscoreweight/sumweight;
Top10scoreweight = Top10scoreweight/sumweight;
XCorrweight = XCorrweight/sumweight;


enscore = Top10score*Top10scoreweight + pscore*Pscoreweight + ...
    perIonmatch*PerIonMatchweight + norm_xcorr*XCorrweight;

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

% glycovweight = .1;
% pepcovweight = .35;
% pscoreweight = .1;
% top10scoreweight = 1/4;
% xcorrweight  = .1;
% Oxocovweight = .1;
Oxoweight = GlycoPATConstant.Scoreweight_EThcD_Oxo;
Y012weight = GlycoPATConstant.Scoreweight_EThcD_Y012;
Glycovweight = GlycoPATConstant.Scoreweight_EThcD_Glycov;
Pepcovweight = GlycoPATConstant.Scoreweight_EThcD_Pepcov;
PerIonMatchweight = GlycoPATConstant.Scoreweight_EThcD_PerIonMatch;
Pscoreweight = GlycoPATConstant.Scoreweight_EThcD_Pscore;
Top10scoreweight = GlycoPATConstant.Scoreweight_EThcD_Top10;
XCorrweight = GlycoPATConstant.Scoreweight_EThcD_XCorr;
sumweight = Oxoweight + Y012weight + Glycovweight + Pepcovweight + PerIonMatchweight + ...
    Pscoreweight + Top10scoreweight + XCorrweight;
Oxoweight = Oxoweight/sumweight;
Y012weight = Y012weight/sumweight;
Glycovweight = Glycovweight/sumweight;
Pepcovweight = Pepcovweight/sumweight;
PerIonMatchweight = PerIonMatchweight/sumweight;
Pscoreweight = Pscoreweight/sumweight;
Top10scoreweight = Top10scoreweight/sumweight;
XCorrweight = XCorrweight/sumweight;


if Top10 >= 5
    Top10score = 1;
else
    Top10score = Top10/5;
end

enscore = glycov*Glycovweight + pepcov*Pepcovweight + ...
    pscore*Pscoreweight + norm_xcorr*XCorrweight + ...
    Top10score*Top10scoreweight + Oxocov*Oxoweight;

if(isnan(enscore))
    enscore =0;
end
if(enscore<0)
    enscore = 0;
end
end