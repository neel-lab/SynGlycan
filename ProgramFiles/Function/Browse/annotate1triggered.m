function annotate1triggered(msdata,handles,scoretable,currentind)
% ANNOTATE1TRIGGERED: annotate 1 set of triggered data (e.g. HCD -> CID,
% ETHCD)
% 
% Syntax:
% annotate1triggered(msdata,handles,scoretable,currentind)
%
% Input:
% msdata: structure, compiled experimental data.
% handles: BrowseGUI main handle, contains information necessary for
%     annotation.
% scoretable: result data, originally a structure, converted to table.
% currentind: which line of the scoretable to be displayed.
%
% Output:
% A parent figure with multiple annotated graph and a table showing scores.
%     Each graph has a button below which generates a more detailed
%     annotation figure.
%
% Note:
% N/A
%
% Example:
% Set breakpoint at the beginning of this program, click 1 cell when
%     browsing Exampleresult_regular.mat with BROWSEGUI
%
% Children function: 
% ANNOTATE1SPECTRUM
%

if length(currentind) > 1
    tgA = figure('Units','Pixels','Position',[50 50 1024 728]);  % triggeredAnnotate / tgA
    spectra = msdata.spectra;
    scannum = msdata.scannum;
    charge = msdata.charge;
    selectedtbl = scoretable(currentind,:);
    numannofig = size(selectedtbl,1);
    uiscoretable = uitable('Parent',tgA,'Unit','normalized','Position',[.05 .8 .9 .17]);
    disptabledata = table2cell(selectedtbl);
    for ii = 1:size(disptabledata,2)
        if isnumeric(disptabledata{1,ii}) && length(disptabledata{1,ii}) > 1
            for j = 1:size(disptabledata,1)
                disptabledata{j,ii} = num2str(disptabledata{j,ii});
            end
        end
    end
    set(uiscoretable,'ColumnName',selectedtbl.Properties.VariableNames,'Data',disptabledata,...
        'columnwidth','auto');
    figposfactor = 0.95/numannofig;
    sgp = selectedtbl.SGP{1};
    theores = cell(numannofig,1);
    theofrag = cell(numannofig,1);
    spectrum = cell(numannofig,1);
    allscannum = selectedtbl.Scan;
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
    thisfigsto = cell(1,numannofig);
    annobuttonsto = cell(1,numannofig);
    
    compid = str2num(selectedtbl.ProteinID{1});
    fragmethodind = ismember(upper(allfragmethods),upper(selectedtbl.Fragmode));
    gettheofragoptions.pepstomode = 0;  % REGENERATE PEP FRAG
    gettheofragoptions.gpfragstomode = false;  % DON'T STORE GP FRAG
    gettheofragoptions.monosacislabile = scoreintdata.sliminput.monosacislabile;
    gettheofragoptions.simultaneousfrag = scoreintdata.sliminput.simultaneousfrag;
    gettheofragoptions.maxstublen = scoreintdata.sliminput.maxstublen;
    gettheofragoptions.sgpseqmode = 1;
    [theofrag,~,~,~] = gettheofragsto([],fragmethodind,{compid},prot,ptmseq,...
        ptmfragsto,ptmfragstublen,ptmtype,ptmmass,[],[],[],allfragmethods,...
        fragnum,userpeptiontyp,gettheofragoptions);
    
    
    localfragmode = selectedtbl.Fragmode;
    [~,ind1] = sort(localfragmode);
    [~,ind2] = sort(allfragmethods);
    [~,ind3] = sort(ind1);
    ind = ind2(ind3);
    ms2tol = scoreintdata.sliminput.ms2tol(ind);
    ms2tolunit = scoreintdata.sliminput.ms2tolunit(ind);
    theofrag = theofrag(ind);
    
    cisoptions.maxlag = scoreintdata.sliminput.maxlag;
    cisoptions.selectpeak = scoreintdata.sliminput.selectpeak;
    
    for ii = 1:numannofig
        thisscannum = selectedtbl.Scan(ii);
        thisfig = axes('Parent',tgA,'Unit','normalized','Position',[0.05+figposfactor*(ii-1),.1,figposfactor*.85,.65]);
        title(thisfig,selectedtbl.Fragmode{ii});
        xlabel(thisfig,'m/z');
        ylabel(thisfig,'Relative Intenisty')
        hold(thisfig,'on');
        tempspec = spectra{scannum == thisscannum};
        [~,~,ind] = unique(tempspec(:,1));
        thisspec = zeros(max(ind),2);
        for j = 1:max(ind)
            tempmz = tempspec(ind == j,1);
            thisspec(j,:) = [tempmz(1),sum(tempspec(ind == j,2))];
        end
        thischarge = charge(scannum == thisscannum);
        fragmethod = upper(selectedtbl.Fragmode{ii});
        tempprecmh = glypepMW(sgp) + 1.007825032;
        if any(strcmpi(fragmethod,{'ETD','ETCID','ETHCD'}))
            thisspec = removePrecursorIon(thisspec,tempprecmh,thischarge,ms2tol(ii),ms2tolunit{ii});
        end
        thisspec = PolishSpectra(thisspec,handles.headerpar.CutOffMed,handles.headerpar.FracMax);
        spectrum{ii} = thisspec;
        
        theores{ii} = calcithscore(thisspec,[theofrag{ii}.mz],thischarge,...
            ms2tol(ii),ms2tolunit{ii},2,cisoptions);
        matchind = ~cellfun(@isempty,theores{ii}.peakmatchindex);
        bar(thisfig,thisspec(~matchind,1),thisspec(~matchind,2)/max(thisspec(:,2))*100,...
            'edgecolor','r','FaceColor','r');
        bar(thisfig,thisspec(matchind,1),thisspec(matchind,2)/max(thisspec(:,2))*100,...
            'edgecolor','g','FaceColor','g');
        thisfigsto{ii} = thisfig;
        annobuttonsto{ii} = uicontrol('Parent',tgA,'Style','pushbutton','String',...
            'Annotate','Unit','normalized',...
            'Position',[.05+figposfactor*(ii-.7) .01 figposfactor/3 .03],...
            'Callback',{@annobutton_callback,ii});
    end
else
    annotate1spectrum(msdata,handles,scoretable,currentind)
end
    function annobutton_callback(src,event,x)
    genannofig(allscannum(x),spectrum{x},table2struct(selectedtbl(x,:)),...
        ms2tol(x),ms2tolunit{x},...
        msdata,theores{x},theofrag{x},'');
    end
end