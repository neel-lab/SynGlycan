function fragments = fragpep_simp(AAseq,npFrag,iontype)
% FRAGPEP_SIMP: a simplified peptide fragmentation process
% 
% Syntax:
% fragments = fragpep_simp(AAseq,npFrag,iontype)
% 
% Input:
% AAseq: amino acid sequence of peptide.
% npFrag: number of peptide bond cleavages allowed
% iontyp: string, should contain 'b', 'c', 'y', 'z' and 'i'. These are the
% peptide ion types allowed in output.
% 
% Output:
% fragments: 1 x n structure, the theoretical fragments.
% 
% Note:
% Only amino acid sequences are accepted as input. No PTM allowed.
%
% Example:
% theofrag = fragpep_simp('LCPDCPLLAPLNDSR',2,'ibcyz')
% 
% theofrag = 
% 
%   1×148 struct array with fields:
% 
% Children function: 
% N/A

% AAseq is peptide amino acid sequence only, nothing more
defaultiontype = 'bcyzi';
numAA = length(AAseq);
numpepbond = numAA - 1;
allstartend = [];
allfragsgps = {AAseq};  % the "none" fragment
allnpFrag = {0};
allngFrag = {0};
allnmFrag = {0};
allfragmz = {pepMW(AAseq,1) + 1.007825032};
allfragtyp = {'none'};
allfragchg = {1};
glyind = zeros(0,1);
modind = zeros(1,0);
allfragunitind = {{1:numAA,glyind,modind}};
allfragstublen = {zeros(1,0)};
todoiontyp = intersect(lower(iontype),defaultiontype);
if ismember('b',lower(iontype)) && ~ismember('y',lower(iontype))
    error('FRAGPEP_SIMP: Input "iontype" is not in pairs.')
end
if ismember('c',lower(iontype)) && ~ismember('z',lower(iontype))
    error('FRAGPEP_SIMP: Input "iontype" is not in pairs.')
end
for i = 1:npFrag
    fragmentationpos = combnk(1:numpepbond,i);
    startend = [];
    for j = 1:size(fragmentationpos,1)
        thisplan = fragmentationpos(j,:);
        startind = 1;
        while ~isempty(thisplan)
            startend = [startend;[startind,thisplan(1)]];
            startind = thisplan(1) + 1;
            thisplan(1) = [];
        end
        startend = [startend;[startind,numpepbond + 1]];
    end
    isduplicate = false(size(startend,1),1);
    for j = 1:size(startend,1)
        for k = 1:size(allstartend,1)
            if isequal(startend(j,:),allstartend(k,:))
                isduplicate(j) = true;
                break
            end
        end
    end
    startend = startend(~isduplicate,:);
    allstartend = [allstartend;startend];
    fragsgps = {};
    fragmzs = {};
    fragtyp = {};
    fragunitind = {};
    fragstublen = {};
    for j = 1:size(startend,1)
        thisfragsgps = AAseq(startend(j,1):startend(j,2));
        % calculate mz
        [~,pepB,pepY,pepC,pepZ] = pepMW(thisfragsgps,1);
        if startend(j,1) == 1
            if ismember('b',todoiontyp)
                fragsgps = [fragsgps,thisfragsgps];
                fragtyp = [fragtyp,['b',num2str(startend(j,2) - startend(j,1) + 1)]];
                fragmzs = [fragmzs,pepB];
                fragunitind = [fragunitind,{{startend(j,1):startend(j,2),glyind,modind}}];
                fragstublen = [fragstublen,{zeros(1,0)}];
            end
            if ismember('c',todoiontyp)
                fragsgps = [fragsgps,thisfragsgps];
                fragtyp = [fragtyp,['c',num2str(startend(j,2) - startend(j,1) + 1)]];
                fragmzs = [fragmzs,pepC];
                fragunitind = [fragunitind,{{startend(j,1):startend(j,2),glyind,modind}}];
                fragstublen = [fragstublen,{zeros(1,0)}];
            end
        elseif startend(j,2) == numAA
            if ismember('y',todoiontyp)
                fragsgps = [fragsgps,thisfragsgps];
                fragtyp = [fragtyp,['y',num2str(startend(j,2) - startend(j,1) + 1)]];
                fragmzs = [fragmzs,pepY];
                fragunitind = [fragunitind,{{startend(j,1):startend(j,2),glyind,modind}}];
                fragstublen = [fragstublen,{zeros(1,0)}];
            end
            if ismember('z',todoiontyp)
                fragsgps = [fragsgps,thisfragsgps];
                fragtyp = [fragtyp,['z',num2str(startend(j,2) - startend(j,1) + 1)]];
                fragmzs = [fragmzs,pepZ];
                fragunitind = [fragunitind,{{startend(j,1):startend(j,2),glyind,modind}}];
                fragstublen = [fragstublen,{zeros(1,0)}];
            end
        else
            if ismember('i',todoiontyp) && ismember('b',todoiontyp) && ismember('c',todoiontyp)
                fragsgps = [fragsgps,thisfragsgps];
                fragtyp = [fragtyp,['i',num2str(startend(j,1)),'-',num2str(startend(j,2))]];
                fragmzs = [fragmzs,pepB];
                fragunitind = [fragunitind,{{startend(j,1):startend(j,2),glyind,modind}}];
                fragstublen = [fragstublen,{zeros(1,0)}];        
                fragsgps = [fragsgps,thisfragsgps];
                
             %% BELOW IS A SEGMENT THAT CONSIDERS IRREGULAR INTERMEDIATE FRAG.
                % HIDDEN FOR CURRENT VERSION
%                 fragtyp = [fragtyp,['iyc',num2str(startend(j,1)),'-',num2str(startend(j,2))]];
%                 fragmzs = [fragmzs,pepB + 17.026549101];
%                 fragunitind = [fragunitind,{{startend(j,1):startend(j,2),glyind,modind}}];
%                 fragstublen = [fragstublen,{zeros(1,0)}];
%                 fragsgps = [fragsgps,thisfragsgps];
%                 fragtyp = [fragtyp,['izb',num2str(startend(j,1)),'-',num2str(startend(j,2))]];
%                 fragmzs = [fragmzs,pepB - 17.026549101];
%                 fragunitind = [fragunitind,{{startend(j,1):startend(j,2),glyind,modind}}];
%                 fragstublen = [fragstublen,{zeros(1,0)}];        
            else
                fragsgps = [fragsgps,thisfragsgps];
                fragtyp = [fragtyp,['i',num2str(startend(j,1)),'-',num2str(startend(j,2))]];
                fragmzs = [fragmzs,pepB];
                fragunitind = [fragunitind,{{startend(j,1):startend(j,2),glyind,modind}}];
                fragstublen = [fragstublen,{zeros(1,0)}];        
            end
        end
    end
    fragnpFrag = num2cell(ones(size(fragmzs)) * i);
    fragngFrag = num2cell(zeros(size(fragmzs)));
    fragnmFrag = num2cell(zeros(size(fragmzs)));
    fragcharge = num2cell(ones(size(fragmzs)));
    allfragsgps = [allfragsgps,fragsgps];
    allnpFrag = [allnpFrag,fragnpFrag];
    allngFrag = [allngFrag,fragngFrag];
    allnmFrag = [allnmFrag,fragnmFrag];
    allfragmz = [allfragmz,fragmzs];
    allfragtyp = [allfragtyp,fragtyp];
    allfragchg = [allfragchg,fragcharge];
    allfragunitind = [allfragunitind,fragunitind];
    allfragstublen = [allfragstublen,fragstublen];
end
% fragments = struct('original',AAseq,'sgp',allfragsgps,'npFrag',allnpFrag,...
%     'ngFrag',allngFrag,'nmFrag',allnmFrag,'mz',allfragmz,...
%     'type',allfragtyp,'charge',allfragchg,'unitindex',allfragunitind,...
%     'stublen',allfragstublen);
fragments = struct('original',AAseq,'sgp',allfragsgps,'npFrag',allnpFrag,...
    'ngFrag',allngFrag,'nmFrag',allnmFrag,'mz',allfragmz,...
    'type',allfragtyp,'charge',allfragchg,'unitindex',allfragunitind);
end