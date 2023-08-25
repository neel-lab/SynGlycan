function annotate1spectrum(msdata,handles,scoretable,currentind)
% ANNOTATE1SPECTRUM: this function is written to provide necessary info for
% GENANNOFIG 
%
% Syntax:
% annotate1spectrum(msdata,handles,scoretable,currentind)
%
% Input:
% msdata: structure, compiled experimental data.
% handles: BrowseGUI main handle, contains information necessary for
%     annotation.
% scoretable: result data, originally a structure, converted to table.
% currentind: which line of the scoretable to be displayed.
% 
% Output:
% Annotation figure
%
% Note:
% N/A
%
% Example:
% Set breakpoint at the beginning of this program, click 1 cell when
%     browsing Exampleresult_regular.mat with BROWSEGUI
%
% Children function: 
% N/A
%

spectra = msdata.spectra;
scannum = msdata.scannum;
charge = msdata.charge;
selectedtbl = scoretable(currentind,:);
sgp = selectedtbl.SGP{1};
scoreintdata = handles.scoreintdata;
prot = scoreintdata.prot;
ptmfragsto = scoreintdata.ptmfragsto;
ptmfragstublen = scoreintdata.ptmfragstublen;
ptmseq = scoreintdata.ptmseq;
ptmtype = scoreintdata.ptmtype;
ptmmass = scoreintdata.ptmmass;
fragnum = scoreintdata.sliminput.fragnum;
userpeptiontyp = scoreintdata.scoreoptions.userpepiontyp;

allfragmethods = scoreintdata.sliminput.fragmode;

compid = str2num(selectedtbl.ProteinID{1});
thisfragmode = selectedtbl.Fragmode{1};
fragmethodind = ismember(upper(allfragmethods),upper(thisfragmode));
gettheofragoptions.pepstomode = 0;  % REGENERATE PEP FRAG
gettheofragoptions.gpfragstomode = false;  % DON'T STORE GP FRAG
gettheofragoptions.monosacislabile = scoreintdata.sliminput.monosacislabile;
gettheofragoptions.simultaneousfrag = scoreintdata.sliminput.simultaneousfrag;
gettheofragoptions.maxstublen = scoreintdata.sliminput.maxstublen;
gettheofragoptions.sgpseqmode = 1;
[theofrag,~,~,~] = gettheofragsto([],fragmethodind,{compid},prot,ptmseq,...
    ptmfragsto,ptmfragstublen,ptmtype,ptmmass,[],[],[],allfragmethods,...
    fragnum,userpeptiontyp,gettheofragoptions);
theofrag = theofrag{1};

ms2tol = scoreintdata.sliminput.ms2tol(fragmethodind);
ms2tolunit = scoreintdata.sliminput.ms2tolunit{fragmethodind};

cisoptions.maxlag = scoreintdata.sliminput.maxlag;
cisoptions.selectpeak = scoreintdata.sliminput.selectpeak;

thisscannum = selectedtbl.Scan(1);
tempspec = spectra{scannum == thisscannum};
[~,~,ind] = unique(tempspec(:,1));
thisspec = zeros(max(ind),2);
for ii = 1:max(ind)
    tempmz = tempspec(ind == ii,1);
    thisspec(ii,:) = [tempmz(1),sum(tempspec(ind == ii,2))];
end
thischarge = charge(scannum == thisscannum);
tempprecmh = glypepMW(sgp) + 1.007825032;
if any(strcmpi(thisfragmode,{'ETD','ETCID','ETHCD'}))
    thisspec = removePrecursorIon(thisspec,tempprecmh,thischarge,ms2tol,ms2tolunit);
end
thisspec = PolishSpectra(thisspec,handles.headerpar.CutOffMed,handles.headerpar.FracMax);
theores = calcithscore(thisspec,[theofrag.mz],thischarge,...
    ms2tol,ms2tolunit,2,cisoptions);
genannofig(thisscannum,thisspec,table2struct(selectedtbl),...
    ms2tol,ms2tolunit,msdata,theores,theofrag,'');
end