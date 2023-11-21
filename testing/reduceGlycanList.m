% Commnets: Program used to reduce glycan search library from 1008 to 928 by removing carbohydrate structure that lack unique 
% diagnostic ions. There are 1006 glycans in the original search library
% that was then reduced to 916.
clc
clear
%  read the text file containing glycans from VPTM file
glycan_txt=readtable('C:\Users\neel\Documents\GitHub\SynGlycan\ProgramFiles\dataFolder\VPTM_N1008.txt')
glycan=glycan_txt(2:end,1)
glycan=table2cell(glycan)
glycan_txt=table2cell(glycan_txt);

% Load allfragments previous generated, with 1008 PTM modifications
allfragments=load('C:\Users\neel\Documents\GitHub\SynGlycan\allfragments_N1008.mat');
allfragments=allfragments.allfragments

% for each glycan calculate in glycan_text file, calculate the precursor
% M+H mass and also all the fragments formed upon fragmentation with ngFrag
% = 1
for i=1:length(glycan)
    frag{i}=multiSGPFrag(glycan{i},0,1,0,1);
    mass{i}=[frag{i}.mz];
    structure{i} = frag{i}(1).original;
    precursor(i) = frag{i}(1).mz;
end

% go through each glycan starting from the last one
for i=length(glycan):-1:1
    idx = find(precursor(1:i-1) == precursor(i));
    % if precursor masses are same for multiple glycans then we may have to
    % delete the glycan provided other conditions are met
    if length(idx)>0
        for j=1:length(idx)
            array_candidate = sort(mass{i});
            array_index = sort(mass{idx(j)});
            % if the number of fragments is same then we may have to
            % deleted the glycan provided the fragment masses are the same
            if length(array_candidate)==length(array_index)
                isEqualOverlap = isequal(array_index, array_candidate(1:length(array_index)));
                % if all fragments are same then isEqualOverlap == 1
                if all(isEqualOverlap==1)
%                    drawglycan(glycan{i},'inputformat','SGP1');
%                    drawglycan(glycan{idx(j)},'inputformat','SGP1');
                % find the index in all fragments that needs to be deleted
                    idx1 = find(strcmp([allfragments.ptmseq], glycan{i}))
                % delete item in glycan list
                    glycan(i)=[];
                    precursor(i)=[];
                    glycan_txt(i+1,:)=[];
                % delete item in allfragments 
                    allfragments.ptmseq(idx1)=[]
                    allfragments.ptmmass(idx1)=[]
                    allfragments.ptmmassfix(idx1)=[]
                    allfragments.ptmtype(idx1)=[]
                    allfragments.ptmfragsto(:,idx1)=[]
                    allfragments.ptmfragstublen(:,idx1)=[]
                    allfragments.ptmisvar(idx1)=[]
                    allfragments.ptmisglycan(idx1)=[]
%                   close all
                    break
                end
            end
        end
    end
end
% save new glycan list as space separated file, independently also save the
% revised allfragments library in .mat format for use in
% GlycoPAT2/SynGlycan
writetable(cell2table(glycan_txt), 'C:\Users\neel\Documents\GitHub\SynGlycan\ProgramFiles\dataFolder\VPTM_N926.txt', 'Delimiter', ' ');
a=1