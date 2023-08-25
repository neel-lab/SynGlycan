function theofrag = fragpep(pepseq,nfrag,iontyp,fragmode)
% FRAGPEP: Perform peptide fragmentation
% 
% Syntax:
% theofrag = fragpep(pepseq,nfrag,iontyp,fragmode)
% 
% Input:
% pepseq: string, amino acid sequence with "{...}" style markers
%     representing PTMs.
% nfrag: 1 x 3 numerical array. [npFrag, ngFrag, nmFrag]. Only npFrag is
%     used.
% iontyp: string, should contain 'b', 'c', 'y', 'z' and 'i'. These are the
%     peptide ion types allowed in output.
% fragmode: string, controls fragmentation behavior when isobaric tags are
%     involved.
% 
% Output:
% theofrag: 1 x n structure, peptide fragments, peptide mass only, all PTMs
%     are treated as zero mass
% 
% Note:
% N/A
%
% Example:
% theofrag = fragpep('LC{1}PDC{1}PLLAPLN{2}DSR',[1,0,0],'bcyz','HCD')
% 
% theofrag = 
% 
%   1×57 struct array with fields:
% 
% theofrag = fragpep('LC{1}PDC{1}PLLAPLN{2}DSR',[2,0,0],'ibcyz','HCD')
% 
% theofrag = 
% 
%   1×148 struct array with fields:
% 
%
% Children function: 
% FRAGPEP_SIMP

% IN DEVELOPMENT - ISOBARIC TAG SUPPORT TO BE ADDED
npfrag = nfrag(1);
[p,g,m] = breakGlyPep(pepseq);
if ~isempty(m)
    allmods = {m.struct};
    isotagind = ismember(allmods,Modification.isotagname);
else
    isotagind = 0;
end
if any(isotagind)  % isobaric tag present (NOT IN USE - RESERVED)
    switch upper(fragmode)
        case {'ETD','ETCID','ETHCD'}
            
        case 'HCD'
            
    end
    
else  % no isobaric tag
    if npfrag == 0
        theofrag.original = pepseq;
        theofrag.sgp = pepseq;
        theofrag.npFrag = 0;
        theofrag.ngFrag = 0;
        theofrag.nmFrag = 0;
        theofrag.mz = pepMW(p.pep,1) + 1.007825032;
        theofrag.type = 'none';
        theofrag.charge = 1;
        theofrag.unitindex = {1:length(p.pep),zeros(0,1),zeros(0,1)};
    elseif npfrag > 0
        theofrag = fragpep_simp(p.pep,npfrag,iontyp);  % backbone mass is calculated here
        % FRAGPEP_SIMP is so simple that it cannot handle any PTM at all,
        % combination with PTM is done below
        ptmstruct = cell(length(g)+length(m),1);
        ptmpos = zeros(length(g)+length(m),1);
        ii = 0;
        if ~isempty(g)
            for ii = 1:length(g)
                ptmstruct{ii} = g(ii).struct;
                ptmpos(ii) = g(ii).pos;
            end
        end
        if ~isempty(m)
            for j = 1:length(m)
                ptmstruct{ii+j} = m(j).struct;
                ptmpos(ii+j) = m(j).pos;
            end
        end
        pepunitindminmax = zeros(length(theofrag),2);  % find what PTM should be on pep frag
        for ii = 1:length(theofrag)
            pepunitindminmax(ii,:) = [min(theofrag(ii).unitindex{1}),max(theofrag(ii).unitindex{1})];
        end
        for ii = 1:length(theofrag)
            onpepptm = ptmpos >= pepunitindminmax(ii,1) & ptmpos <= pepunitindminmax(ii,2);
            if any(onpepptm)
                fragpepseq = theofrag(ii).sgp;
                onpepptmstruct = ptmstruct(onpepptm);
                onpepptmpos = ptmpos(onpepptm) - pepunitindminmax(ii,1) + 1;
                for j = length(onpepptmpos):-1:1
                    fragpepseq = [fragpepseq(1:onpepptmpos(j)),onpepptmstruct{j},fragpepseq(onpepptmpos(j)+1:end)];
                end
                theofrag(ii).sgp = fragpepseq;  % This is the final pep frag - PTM mass is ignored
            end
            if ii == 1
                [theofrag.original] = deal(theofrag(ii).sgp);  % because 1st frag ion is always "none"
            end
        end
    end
end
end