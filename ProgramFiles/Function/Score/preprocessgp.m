function [ptminfo,prot] = preprocessgp(prot,varptm)
% PREPROCESSGP: read digested SGP file, convert fixed PTM to variable PTM
% like ones
% 
% Syntax:
% [ptminfo,prot] = preprocessgp(prot,varptm)
% 
% Input:
% prot: m x 1 cell, each element is a n x 1 cell array of strings. Each m
% is a protein, each n is a (glyco)peptide
% 
% Output:
% ptminfo: structure with field "mod" and "uniind". "mod" is PTM structure,
% "uniind" is their indices.
% prot: same structure as input "prot", but all glycopeptides were
% converted.
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

glycan1let = Glycan.gly1let;
stoptm_struct = varptm.mod;
stoptm_ind = varptm.varptmuniind;
newptm_startind = max(stoptm_ind) + 1;
for i = 1:length(prot)
    % i means each protein
    thisprot = prot{i};
    %% diagnose
    %     thisprot_ref = thisprot;
    %% diagnose end
    for j = 1:length(thisprot)
        % j means each (glyco)peptide
        thispep = thisprot{j};
        [p,g,m] = breakGlyPep(thispep);
        % p: pure peptide AA
        % g: variable PTM, or glycan as fixed PTM
        % m: fixed non-glycan PTM only
        % separate glycans inside g
        
        %% VERSION 1 - ORIGINAL DESIGN
        ptmtype = [];
        % 1 for non-glycan PTM, fixed
        % 2 for non-glycan PTM, variable
        % 3 for glycan PTM, fixed
        % 4 for glycan PTM, variable
        ptmpos = [];
        ptmstruct = {};
        if ~isempty(m)  % fixed non-glycan PTM
            for k = 1:length(m)
                ptmtype = [ptmtype;1];
                ptmpos = [ptmpos;m(k).pos];
                ptmstruct = [ptmstruct;m(k).struct];
            end
        end
        if ~isempty(g)  % fixed or variable glycan PTM, variable non-glycan PTM
            for k = 1:length(g)
                % determine if it's {+/-XXX.XX} or {X}, X marks a digit
                if ismember(g(k).struct,{'+','-'})  % fixed glycan PTM, represented by its mass
                    ptmtype = [ptmtype;3];
                    ptmstruct = [ptmstruct;g(k).struct];
                elseif ismember(g(k).struct,glycan1let)  % fixed glycan PTM, represented by SGP1
                    ptmtype = [ptmtype;3];
                    ptmstruct = [ptmstruct;g(k).struct];
                else  % below are variable PTMs specified in digest file
                    thisvarptmind = str2double(g(k).struct(2:end-1));
                    thesemod = stoptm_struct(stoptm_ind == thisvarptmind);
                    if strcmp(thesemod{1}(1),'{')  % vartpm: glycan
                        ptmtype = [ptmtype;4];
                        ptmstruct = [ptmstruct;g(k).struct];
                    elseif strcmp(thesemod{1}(1),'<')  % varptm: non-glycan
                        ptmtype = [ptmtype;2];
                        ptmstruct = [ptmstruct;g(k).struct];
                    end
                end
                ptmpos = [ptmpos;g(k).pos];
            end
        end
        %% VERSION 1 END
        
        %% VERSION 2
%         ptmtype = zeros(1,length(g) + length(m));
%         ptmpos = zeros(1,length(g) + length(m));
%         ptmstruct = cell(1,length(g) + length(m));
%         ptmcounter = 0;
%         if ~isempty(m)  % fixed non-glycan PTM
%             for k = 1:length(m)
%                 ptmtype(k) =1;
%                 ptmpos(k) = m(k).pos;
%                 ptmstruct{k} = m(k).struct;
%             end
%             ptmcounter = ptmcounter + k;
%         end
%         if ~isempty(g)  % fixed or variable glycan PTM, variable non-glycan PTM
%             for k = 1:length(g)
%                 % determine if it's {+/-XXX.XX} or {X}, X marks a digit
%                 if ismember(g(k).struct,{'+','-'})  % fixed glycan PTM, represented by its mass
%                     ptmtype(k + ptmcounter) = 3;
%                     ptmstruct{k + ptmcounter} = g(k).struct;
%                 elseif ismember(g(k).struct,glycan1let)  % fixed glycan PTM, represented by SGP1
%                     ptmtype(k + ptmcounter) = 3;
%                     ptmstruct{k + ptmcounter} = g(k).struct;
%                 else  % below are variable PTMs specified in digest file
%                     thisvarptmind = str2double(g(k).struct(2:end-1));
%                     thesemod = stoptm_struct(stoptm_ind == thisvarptmind);
%                     if strcmp(thesemod{1}(1),'{')  % vartpm: glycan
%                         ptmtype(k + ptmcounter) = 4;
%                         ptmstruct = g(k).struct;
%                     elseif strcmp(thesemod{1}(1),'<')  % varptm: non-glycan
%                         ptmtype(k + ptmcounter) = 2;
%                         ptmstruct{k + ptmcounter} = g(k).struct;
%                     end
%                 end
%                 ptmpos(k + ptmcounter) = g(k).pos;
%             end
%         end
        %% VERSION 2 END
        
        isvarptm = (ptmtype == 2) | (ptmtype == 4);
        fixptmstruct = ptmstruct(~isvarptm);
        fixptmpos = ptmpos(~isvarptm);
        ispresent = ismember(fixptmstruct,stoptm_struct);
        exist_ptmstruct = fixptmstruct(ispresent);
        exist_ptmpos = fixptmpos(ispresent);
        exist_ptmnewind = zeros(size(exist_ptmpos));
        for k = 1:length(exist_ptmstruct)
            exist_ptmnewind(k) = stoptm_ind(ismember(stoptm_struct,exist_ptmstruct{k}));
        end
        add_ptmstruct = ptmstruct(~ispresent);
        add_ptmpos = ptmpos(~ispresent);
        [uniadd_ptmstruct,~,ind] = unique(add_ptmstruct);
        add_ptmnewind = (newptm_startind:newptm_startind + max(ind) - 1)';
        if any(~ispresent)
            stoptm_struct = [stoptm_struct;uniadd_ptmstruct];
            stoptm_ind = [stoptm_ind;add_ptmnewind];
        end
        for k = 1:length(exist_ptmstruct)
            exist_ptmstruct{k} = ['{',num2str(exist_ptmnewind(k)),'}'];
        end
        for k = 1:max(ind)
            [add_ptmstruct{ind == k}] = deal(['{',num2str(add_ptmnewind(k)),'}']);
        end
        newptm = struct('struct',[ptmstruct(isvarptm);exist_ptmstruct;add_ptmstruct],'pos',num2cell([ptmpos(isvarptm);exist_ptmpos;add_ptmpos]));
        [~,ind]=sort([newptm.pos]);
        newptm = newptm(ind);
        thisprot{j} = joinGlyPep(p,newptm,[]);
        
        newptm_startind = newptm_startind + sum(~ispresent);
    end
    %% diagnose
    %     prot_comp = thisprot;
    %     for j = 1:length(thisprot)
    %         for k = max(stoptm_ind) + 1:length(stoptm_ind)
    %             thisprot{j} = strrep(thisprot{j},['{',num2str(k),'}'],stoptm_struct{k});
    %         end
    %     end
    %     if ~isequal(prot_comp,thisprot)
    %         i
    %     end
    %% diagnose end
    prot{i} = thisprot;
end

ptminfo.mod = stoptm_struct;
ptminfo.uniind = stoptm_ind;
end