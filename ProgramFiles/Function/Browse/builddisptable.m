function [displaydata,displaydataind,displaycolnames,displaycolformat] = ...
    builddisptable(scoredata,displaydataind,displayoptions,usercustom)
% BUILDDISPTABLE: create the table for display in BROWSEGUI, the column
%     names of the table and where each identification in the result is
%     placed.
%
% Syntax:
% [displaydata,displaydataind,displaycolnames,displaycolformat] = builddisptable(...
%     scoredata,displayoptions,usercustom)
%
% Input:
% scoredata: 1 x n structure. Scoring results.
% displayoptions: structure. User setting on how to organize the results.
% usercostom: structure. User added information. Currently there are 2
%     fields: a n x 3 logical array "Quality3" and a n x 1 cell array of strings
%     "SGPwas". "Quality3" keeps  record of user's decision on spectrum
%     quality, precursor ion selection and user's subjective decision of
%     every identification. "SGPwas" keeps the original SGP of
%     candidates in case user decides that another candidate fits the
%     spectrum better.
%
% Output:
% displaydata: m x n cell array. The data to be displayed in BROWSEGUI
%     window.
% displaydataind: m x n cell array or double. The serial number of results
%     that's in each cell of "displaydata".
% displaycolnames: 1 x n cell array of strings. The column names of the
%     table. If the experiment is triggered, "displaycolnames" comes in 2
%     parts: basic information and fragmentation method specific ones.
%     Otherwise, the columns will be consisted of 3 criteria and all the
%     parameters in the result.
% displaycolformat: 1 x n cell array of strings. The format of each column.
%     The values can be fount at MATLAB help document of UITABLE.
%
%     (Below are the columns for triggered experiments, depending on
%     whether the isomers are combined format will be different - combined
%     first.)
%     Basic information and formats are:
%     Spectrum quality - logical
%     Monoisotopic precursor ion assignment - logical
%     User subjective decision - logical
%     Protein name - char
%     Precursor ion monoisotopic mass - numeric
%     Candidate theoretical mass - numeric
%     Precursor ion charge - numeric
%     Candidate SGP sequence - char
%     Combined enscore - numeric
%     Combined decoy enscore - numeric
%     AUC - numeric
%     Retention time (precursor ion) - numeric
%     Original candidate SGP sequence - char
%
%     Fragmentation method specific information and formats are:
%     Scan number - char / numeric
%     Peak lag - char / numeric
%     P-value - char / numeric
%     Top 10 - char / numeric
%     % Ion match - char / numeric
%     Enscore - char / numeric
%     Decoy enscore - char / numeric
%
% Note:
% If the experiment is triggered and user chose to combine isomers, the
%     displayed data will be sorted according to combiES.
%
% Example:
% N/A
%
% Children function:
%

% GlycoPAT 2 authors: Kai Cheng, Gabrielle Pawlowski, Sriram Neelamegham
%(c) 2020, Research Foundation for State University of New York. All rights reserved

combineisomer = displayoptions.combineisomer;
scoreintdata = displayoptions.scoreintdata;

protein = {scoredata.Protein};
scans = [scoredata.Scan];
sgps = {scoredata.SGP};
expt = [scoredata.Expt];
theo = [scoredata.Theo];
charge = [scoredata.Charge];
peaklag = [scoredata.PeakLag];
pvalue = [scoredata.Pvalue];
top10 = [scoredata.Top10];
percentionmatch = [scoredata.PercentIonMatch];
enscore = [scoredata.Enscore];
decoyes = [scoredata.DecoyEnscore];
if isfield(scoredata,'localization')
    locals = {scoredata.localization};
else
    locals = cell(size(protein));
end
quant = [scoredata.Quant];
protids = {scoredata.ProteinID};
quality3 = usercustom.Quality3;
SGPwas = usercustom.SGPwas';    %Added {}, removed {} 7/31/2023 TH
%{
if size(SGPwas) == [1,1]
    SGPwas(1,1) = SGPwas;
    SGPwas = SGPwas';
end
%}

if scoreintdata.scoreoptions.isHCDtrigger  % Trigger experiment
    SASSO = scoreintdata.SASSO;
    SASSO_method = scoreintdata.colnames;
    nummethods = length(SASSO_method);
    displaycolnames = {'SpecQual','MS1Qual','Subjective',...
        'Protein','Expt','Theo','Chg','SGP',...
        'SGPwas','CombiES','decoyCombiES','AUC','Precursorscan',...
        'Retime'};
    displaycolformat = {'logical','logical','logical',...
        'char','numeric','numeric','numeric','char',...
        'char','numeric','numeric','numeric','numeric',...
        'numeric'};
    % 9 basic columns
    for ii = 1:nummethods
        tempmethod = SASSO_method{ii};
        displaycolnames = [displaycolnames,['Scan','_',tempmethod],['peakLag','_',tempmethod],...
            ['P','_',tempmethod],['Top10','_',tempmethod],['%IonMatch','_',tempmethod],...
            ['ES','_',tempmethod],['decoyES','_',tempmethod]];  % Each fragmode has 7 columns
    end
    displaycolnames = [displaycolnames,'localization','Comment'];
    if isempty(displaydataind)
        if combineisomer  % Trigger experiment - combine isomer
            displaydataind = grouptriggeredresult(SASSO,protids,scans,'combineisomer');
        else  % Trigger experiment - do not combine isomer
            displaydataind = grouptriggeredresult(SASSO,protids,scans,'regular');
        end
    end
    disptbht = size(displaydataind,1);
    disptable_basic = cell(disptbht,14);  % Default info
    disptable_expand = cell(disptbht,7*nummethods + 1);  % Method specific table
    displaycolformat_expand = cell(1,7*nummethods + 1);
    [CombiES,CombiDecoyES] = calculateCombiES(scoredata,displaydataind,...
        SASSO_method);
    if combineisomer
        for ii = 1:disptbht  % "newtableind" contains zeroes in case scan got no hit
            % i means each isomer group
            tempdispdatainds = displaydataind(ii,:);
            numisomers = length(tempdispdatainds{1});
            isomergrpdispdatainds = zeros(numisomers,nummethods);
            tempdefaultdispdataind = zeros(2,numisomers);
            % Picking a representative for the "default" part of the table.
            % It is the first isomer of each group
            % 1st row is the serial number, 2nd row tells which
            % fragmentation mode
            for jj = 1:numisomers  % each isomer
                for kk = 1:nummethods  % each frag method
                    isomergrpdispdatainds(jj,kk) = tempdispdatainds{kk}(jj);
                    % Each row is a isomer, each column is a frag method
                end
                % Keep a record on combiES
                thisisomerdispdatainds = isomergrpdispdatainds(jj,:);
                % Isomers may not go though all fragmentation modes, find
                % at least one for extracting info about its MS2 scan.
                oneisomerind = find(thisisomerdispdatainds,1);
                % Which frag method has a candidate identified
                tempdefaultdispdataind(1,jj) = thisisomerdispdatainds(oneisomerind);
                tempdefaultdispdataind(2,jj) = oneisomerind;
            end
            % In group sorting -  based on combiES - highest first
            [tempcombiES,ind] = sort(CombiES{ii},'descend');
            tempcombidecoyES = CombiDecoyES{ii}(ind);
            isomergrpdispdatainds = isomergrpdispdatainds(ind,:);
            tempprot = '';
            temptheo = '';
            tempsgp = '';
            tempsgpwas = '';
            templocal = '';
            for jj = 1:numisomers
                % i means each isomer group
                % j means each isomer
                tempprot = [tempprot,protein{tempdefaultdispdataind(1,jj)},' '];
                temptheo = [temptheo,num2str(theo(tempdefaultdispdataind(1,jj))),' '];
                tempsgp = [tempsgp,sgps{tempdefaultdispdataind(1,jj)},' '];
                templocal = [templocal,locals{tempdefaultdispdataind(1,jj)},' '];
                if strcmpi(SGPwas{tempdefaultdispdataind(1,jj)},sgps{tempdefaultdispdataind(1,jj)})  % No change
                    tempsgpwas = [tempsgpwas,'Original',' '];
                else
                    tempsgpwas = [tempsgpwas,SGPwas{tempdefaultdispdataind(1,jj)},' '];
                end
            end
            tempprot = tempprot(1:end-1);
            temptheo = temptheo(1:end-1);
            tempsgp = tempsgp(1:end-1);
            tempsgpwas = tempsgpwas(1:end-1);
            templocal = templocal(1:end-1);
            tempprecscan = SASSO(SASSO(:,tempdefaultdispdataind(2,1)) == ...
                scans(tempdefaultdispdataind(1,1)),end - 1);
            tempretime = SASSO(SASSO(:,tempdefaultdispdataind(2,1)) == ...
                scans(tempdefaultdispdataind(1,1)),end);
            disptable_basic(ii,:) = {true,true,false,...  % SpecQual MS1Qual Subjective
                tempprot,expt(tempdefaultdispdataind(1)),temptheo,...  % Protein Expt Theo
                charge(tempdefaultdispdataind(1)),tempsgp,tempsgpwas,...  % Chg SGP SGPwas
                num2str(tempcombiES(:)'),num2str(tempcombidecoyES(:)'),...  % CombiES decoyCombiES
                quant(tempdefaultdispdataind(1)),tempprecscan(1),...  % AUC Precursorscan
                tempretime(1)};  % AUC Retime
            % Note: "tempprecscan" and "temppretime" has "(1)" because
            %     preprocessing may assign different triggered scan (CID,
            %     ETHCD, etc.) to same triggerind scan (HCD), in which case
            %     there will be 2 precursor scan numbers, but there should
            %     only be 1. The "(1)" is OK as a patch.
            tempquality3 = {true(numisomers,3),true(numisomers,nummethods),...
                false(numisomers,nummethods)};  % SpecQual, MS1Qual: default true; Subjective: default false
            for jj = 1:nummethods
                tempscans = '';
                temppklags = '';
                temppvalues = '';
                temptop10s = '';
                temppercents = '';
                tempenscores = '';
                tempdecoyes = '';
                for kk = 1:numisomers
                    tempind = isomergrpdispdatainds(kk,jj);
                    if tempind > 0
                        tempscans = [tempscans,' ',num2str(scans(tempind))];
                        temppklags = [temppklags,' ',num2str(peaklag(tempind))];
                        temppvalues = [temppvalues,' ',num2str(pvalue(tempind))];
                        temptop10s = [temptop10s,' ',num2str(top10(tempind))];
                        temppercents = [temppercents,' ',num2str(percentionmatch(tempind))];
                        tempenscores = [tempenscores,' ',num2str(enscore(tempind))];
                        tempdecoyes = [tempdecoyes,' ',num2str(decoyes(tempind))];
                        tempquality3{1}(kk,jj) = quality3(tempind,1);
                        tempquality3{2}(kk,jj) = quality3(tempind,2);
                        tempquality3{3}(kk,jj) = quality3(tempind,3);
                    else
                        tempscans = [tempscans,' N/A'];
                        temppklags = [temppklags,' N/A'];
                        temppvalues = [temppvalues,' N/A'];
                        temptop10s = [temptop10s,' N/A'];
                        temppercents = [temppercents,' N/A'];
                        tempenscores = [tempenscores,' N/A'];
                        tempdecoyes = [tempdecoyes,' N/A'];
                    end
                end
                disptable_expand{ii,(jj-1)*7+1} = tempscans;
                disptable_expand{ii,(jj-1)*7+2} = temppklags;
                disptable_expand{ii,(jj-1)*7+3} = temppvalues;
                disptable_expand{ii,(jj-1)*7+4} = temptop10s;
                disptable_expand{ii,(jj-1)*7+5} = temppercents;
                disptable_expand{ii,(jj-1)*7+6} = tempenscores;
                disptable_expand{ii,(jj-1)*7+7} = tempdecoyes;
                displaydataind{ii,jj} = isomergrpdispdatainds(:,jj);
            end
            disptable_expand{ii,jj*7+1} = templocal;
            disptable_expand{ii,jj*7+2} = '';
            % Rules for displaying Quality3: 1. spectrum quality: any FALSE
            %     makes checkbox false; 2. MS1 quality: same as spectrum
            %     quality; 3. user subjective: any TRUE makes checkbox true.
            if any(any(~tempquality3{1}))
                disptable_basic{ii,1} = false;
            end
            if any(any(~tempquality3{2}))
                disptable_basic{ii,2} = false;
            end
            if any(any(tempquality3{3}))
                disptable_basic{ii,3} = true;
            end
        end
        displaydata = [disptable_basic,disptable_expand];
        [displaycolformat_expand{:}] = deal('char');
        displaycolformat = [displaycolformat,displaycolformat_expand];
    else  % Each isomer will be displayed in a new row, this is
        % actually a simplified version of the section above
        for ii = 1:disptbht
            tempdispdatainds = displaydataind(ii,:);
            reftempdispdataind_ind = find(tempdispdatainds,1);
            reftempdispdataind = tempdispdatainds(reftempdispdataind_ind);
            temprecscan = SASSO(SASSO(:,reftempdispdataind_ind) == ...
                scans(reftempdispdataind),end - 1);
            tempretime = SASSO(SASSO(:,reftempdispdataind_ind) == ...
                scans(reftempdispdataind),end);
            disptable_basic(ii,:) = {true,true,true,...
                protein{reftempdispdataind},expt(reftempdispdataind),...
                theo(reftempdispdataind),charge(reftempdispdataind),sgps{reftempdispdataind},...
                SGPwas{reftempdispdataind},CombiES(ii),CombiDecoyES(ii),...
                quant(reftempdispdataind),temprecscan(1),tempretime(1)};
            tempquality3 = {true(1,nummethods),true(1,nummethods),false(1,nummethods)};
            % SpecQual, MS1Qual: default true; Subjective: default false
            for jj = 1:nummethods
                tempind = tempdispdatainds(jj);
                if tempind > 0
                    tempscans = scans(tempind);
                    temppklags = peaklag(tempind);
                    temppvalues = pvalue(tempind);
                    temptop10s = top10(tempind);
                    temppercents = percentionmatch(tempind);
                    tempenscores = enscore(tempind);
                    tempdecoyes = decoyes(tempind);
                    templocal = locals{tempind};
                    if ~isempty(quality3)
                        tempquality3{1}(jj) = quality3(tempind,1);
                        tempquality3{2}(jj) = quality3(tempind,2);
                        tempquality3{3}(jj) = quality3(tempind,3);
                    end
                else
                    tempscans = -1;
                    temppklags = -50;
                    temppvalues = -1;
                    temptop10s = -1;
                    temppercents = -1;
                    tempenscores = -1;
                    tempdecoyes = -1;
                    templocal = '';
                end
                disptable_expand{ii,(jj-1)*7+1} = tempscans;
                disptable_expand{ii,(jj-1)*7+2} = temppklags;
                disptable_expand{ii,(jj-1)*7+3} = temppvalues;
                disptable_expand{ii,(jj-1)*7+4} = temptop10s;
                disptable_expand{ii,(jj-1)*7+5} = temppercents;
                disptable_expand{ii,(jj-1)*7+6} = tempenscores;
                disptable_expand{ii,(jj-1)*7+7} = tempdecoyes;
            end
            disptable_expand{ii,jj*7+1} = templocal;
            disptable_expand{ii,jj*7+2} = '';
        end
        displaydata = [disptable_basic,disptable_expand];
        [displaycolformat_expand{:}] = deal('numeric');
        displaycolformat_expand{end} = 'char';
        displaycolformat = [displaycolformat,displaycolformat_expand];
    end
else  % Regular
    displaycolnames = fieldnames(scoredata);
    displaycolnames = ['SpecQual','MS1Qual','Subjective','SGPwas',displaycolnames(:)','localization','Comment'];
    displaycolformat = cell(size(displaycolnames));
    [displaycolformat{:}] = deal('char');
    displaycolformat(1:3) = {'logical','logical','logical'};
    numericfld = ResultItem.itemisnumeric;
    numcellfld = ResultItem.itemisnumcell;
    stringfld = ResultItem.itemisstring;
    if combineisomer  % Regular experiment - combine isomer
        uniscans = unique(scans);
        synsasso = zeros(length(uniscans),3);
        synsasso(:,1) = uniscans;  % Synthetic SASSO - cheat the program
        if isempty(displaydataind)
            displaydataind = grouptriggeredresult(synsasso,protids,scans,'combineisomer');
        end
        displaydata = cell(length(displaydataind),length(displaycolnames));
        for ii = 1:length(displaydataind)
            displaydata(ii,1:3) = {true,true,true};
            tempdisplaydataind = displaydataind{ii};
            tempenscore = enscore(tempdisplaydataind);
            [~,ind] = sort(tempenscore,'descend');
            tempdisplaydataind = tempdisplaydataind(ind);
            displaydataind{ii} = tempdisplaydataind;
            for jj = 4:length(displaycolnames)
                tempdisplaycolnames = displaycolnames{jj};
                tempdisplaydata = '';
                if ismember(tempdisplaycolnames,numericfld)
                    for kk = 1:length(tempdisplaydataind)
                        tempdisplaydata = [tempdisplaydata,' ;',...
                            num2str(scoredata(tempdisplaydataind(kk)).(tempdisplaycolnames))];
                    end
                    tempdisplaydata = tempdisplaydata(3:end);
                elseif ismember(tempdisplaycolnames,[numcellfld,stringfld])
                    for kk = 1:length(tempdisplaydataind)
                        tempdisplaydata = [tempdisplaydata,' ;',...
                            scoredata(tempdisplaydataind(kk)).(tempdisplaycolnames)];
                    end
                    tempdisplaydata = tempdisplaydata(3:end);
                elseif strcmpi(tempdisplaycolnames,'SGPwas')
                    for kk = 1:length(tempdisplaydataind)
                        tempdisplaydata = [tempdisplaydata,' ;',...
                            SGPwas{tempdisplaydataind(kk)}];
                    end
                    tempdisplaydata = tempdisplaydata(3:end);
                else
                    tempdisplaydata = '';
                end
                displaydata{ii,jj} = tempdisplaydata;
            end
            if ~isempty(quality3)
                tempquality3 = quality3(tempdisplaydataind,:);
                if ~any(any(tempquality3(:,1)))
                    displaydata{ii,1} = false;
                end
                if ~any(any(tempquality3(:,2)))
                    displaydata{ii,2} = false;
                end
                if any(any(tempquality3(:,3)))
                    displaydata{ii,3} = true;
                end
            end
        end
    else  % Regular experiment - do not combine isomer - display original data
        if isempty(displaydataind)
            displaydataind = (1:length(scoredata))';
        end
        tempscoredata = scoredata(displaydataind);
        tempsgpwas = SGPwas(displaydataind);
        tempquality3 = quality3(displaydataind,:);
        displaydata = cell(length(tempscoredata),length(displaycolnames));
        scoredatafieldnames = fieldnames(tempscoredata);
        for ii = 5:length(displaycolnames)  % the first 3 columns are the 3 user decisions
            if ismember(displaycolnames{ii},scoredatafieldnames)
                tempdata = {tempscoredata.(displaycolnames{ii})};
                displaydata(:,ii) = tempdata(:);
                if ismember(displaycolnames{ii},numericfld)
                    displaycolformat{ii} = 'numeric';
                elseif ismember(displaycolnames{ii},numcellfld)
                    displaycolformat{ii} = 'char';
                elseif ismember(displaycolnames{ii},stringfld)
                    displaycolformat{ii} = 'char';
                end
            else
                displaydata(:,ii) = repmat({''},size(displaydata,1),1);
            end
        end
        for ii = 1:size(displaydata,1)
            if strcmpi(tempscoredata(ii).SGP,tempsgpwas{ii})
                displaydata{ii,4} = 'Original';
            else
                displaydata{ii,4} = tempsgpwas{ii};
            end
        end
        if ~isempty(quality3)
            displaydata(:,1:3) = num2cell(tempquality3);
        end
    end
end
end