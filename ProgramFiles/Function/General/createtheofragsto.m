function alltheofrag = createtheofragsto(sgpcomp,fragmodes,scoreintdata)
% CREATETHEOFRAGSTO: create theoretical fragments for input glycopeptides.
%     This funciton is dedicated to browsers and single spectrum
%     annotaiton, as high-throughput scoring requires higher efficiency
%     code.
%
% Syntax:
% alltheofrag = createtheofragsto(sgpcomp,fragmodes,scoreintdata)
%
% Input:
% sgpcomp: n x 1 cell array of 1 x m numbers. n equals to number of
%     glycopeptides. The format for each 1 x m is: ["Protein #", "Peptide #",
%     "PTM 1", "PTM 2",...]. These serial numbers represents the position
%     of corresponding elements in "scoreintdata".
% fragmodes: n x 1 cell array of strings. The fragmentation modes to use
%     when generating theoretical fragments.
% scoreintdata: structure. A storage of intermediate variables for scoring.
%     This structure must contain the following fields: "ptmfragsto",
%     "ptmfragstublen", "ptmseq", "ptmtype", "ptmmass", "prot",
%     "fragmode", "fragiontyp", "fragnum", "maxstublen", "monosacislabile",
%     "simultaneousfrag" and "addoxoniumion". The detail of  each field can
%     be found in function "scoreAllSpectra" and "theoptmfrag".
%
% Output:
% alltheofrag: m x n cell array of 1 x k structure. Theoretical fragments
%     of input glycopeptides. m equals to the number of glycopeptides, n
%     equals to the number of fragmentation modes.
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
% N/A

% GlycoPAT 2 authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2020, Research Foundation for State University of New York. All rights reserved

ptmfragsto = scoreintdata.ptmfragsto;
ptmfragstublen = scoreintdata.ptmfragstublen;
ptmseq = scoreintdata.ptmseq;
ptmtype = scoreintdata.ptmtype;
ptmmass = scoreintdata.ptmmass;
prot = scoreintdata.prot;
allfragmethods = scoreintdata.sliminput.fragmode;
fragiontyp = scoreintdata.sliminput.fragiontyp;
fragnum = scoreintdata.sliminput.fragnum;
maxstublen = scoreintdata.sliminput.maxstublen;
monosacislabile = scoreintdata.sliminput.monosacislabile;
simultaneousfrag = scoreintdata.sliminput.simultaneousfrag;
addoxoniumion = scoreintdata.sliminput.addoxoniumion;

defaultfragiontyp = 'bcyzi';

alltheofrag = cell(length(sgpcomp),length(fragmodes));

for ii = 1:length(sgpcomp)
    thissgpcomp = sgpcomp{ii};
    if ischar(thissgpcomp)
        thissgpcomp = str2num(thissgpcomp);
    end
    pepseq = prot{thissgpcomp(1)}{thissgpcomp(2)};
    for jj = 1:length(fragmodes)
        tempfragmethod = upper(fragmodes{jj});
        fragmethodind = ismember(upper(allfragmethods),tempfragmethod);
        userpepiontyp = defaultfragiontyp(fragiontyp(fragmethodind,:));
        nfrag = fragnum(fragmethodind,:);
        tempnpfrag = nfrag(1);
        ptmfragrownum = nfrag(2);
        if ptmfragrownum == 0
            ptmfragrownum = 1;
        end
        if strcmpi(tempfragmethod,'CID')
            tempPTMseqs = ptmseq(thissgpcomp(3:end));
            totalmonosac = sum(cellfun(@(x) length(strfind(x,'{')),tempPTMseqs));
            if totalmonosac < 6
                tempnpfrag = max(1,nfrag(1));
            end
        end
        pepfrag = fragpep(pepseq,tempnpfrag,userpepiontyp,tempfragmethod);
        if length(thissgpcomp) > 2
            fragopt.stublen = maxstublen;
            fragopt.mode = 1;
            fragopt.fragmode = tempfragmethod;
            fragopt.monosacislabile = monosacislabile(fragmethodind);
            fragopt.simultaneousfrag = simultaneousfrag(fragmethodind);
            fragopt.addoxoniumion = addoxoniumion(fragmethodind);
            alltheofrag{ii,jj} = combitheofrag(pepfrag,ptmfragsto(ptmfragrownum,:),...
                ptmfragstublen(ptmfragrownum,:),pepseq,ptmseq,...
                ptmtype,ptmmass,thissgpcomp,[tempnpfrag,nfrag(2),nfrag(3)],fragopt);
        else
            alltheofrag{ii,jj} = pepfrag;
        end
    end
end
end