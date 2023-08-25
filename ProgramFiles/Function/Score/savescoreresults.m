function savescoreresults(result,scoreoptions,input,PTM,allpepdata,triggerdata,...
    outputdir,outputfilename,protindse,mode)
% SAVESCORERESULTS: compile file names and save score results 
% 
% Syntax:
% savescoreresults(result,scoreoptions,input,PTM,allpepdata,triggerdata,...
%     outputdir,outputfilename,protindse,mode)
% 
% Input:
% result: 
% scoreoptions
% input
% PTM
% allpepdata
% triggerdata
% outputdir
% outputfilename
% protindse
% mode

% 
% Output:
% 
% 
% Note:
% 
% 
% Example:
% 
% 
% Children function: 
% 
% 
% See Also:
% 

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 04/01/2020

maxresultsize = 500000;
switch mode
    case 'save'
        prot = allpepdata.glypep;
        protname = allpepdata.FASTAhead;
        sliminput = input;
        fldnms = fieldnames(sliminput);
        if ismember('handles',fldnms)
            sliminput = rmfield(sliminput,'handles');
        end
        if ismember('allfragments',fldnms)
            sliminput = rmfield(sliminput,'allfragments');
        end
        if ismember('allfragments_presearch',fldnms)
            sliminput = rmfield(sliminput,'allfragments_presearch');
        end
        scoreintdata.scoreoptions = scoreoptions;
        scoreintdata.sliminput = sliminput;
        scoreintdata.ptmisomers = PTM.isomers;
        scoreintdata.prot = prot;
        scoreintdata.protname = protname;
        scoreintdata.filetype = 'MainSearchSave';
%         if isfield('allfragments',input)
            allptmseq = input.allfragments.ptmseq;
            for ii = 1:length(result)
                tempprotid = str2num(result(ii).ProteinID);
                tempsgp = result(ii).SGP;
                [~,g,m] = breakGlyPep(tempsgp);
                tempptmseq = {};
                tempptmpos = [];
                if ~isempty(g)
                    tempptmseq = [tempptmseq,g.struct];
                    tempptmpos = [tempptmpos,g.pos];
                end
                if ~isempty(m)
                    tempptmseq = [tempptmseq,m.struct];
                    tempptmpos = [tempptmpos,m.pos];
                end
                [~,ind] = sort(tempptmpos);
                tempptmseq = tempptmseq(ind);
                for jj = 1:length(tempptmseq)
                    tempprotid(jj + 2) = find(strcmpi(allptmseq,tempptmseq{jj}));
                end
                result(ii).ProteinID = num2str(tempprotid);
            end
%         end
        if scoreoptions.isHCDtrigger
            SASSO = triggerdata.SASSO;
            colnames = triggerdata.colnames;
            colselected = ismember(upper(colnames),upper(scoreoptions.analyzefragmode));
            colnames = colnames(colselected);
            SASSO = [SASSO(:,colselected),SASSO(:,end-1:end)];
            scoreintdata.SASSO = SASSO;
            scoreintdata.colnames = colnames;
        end
        [p,f,e] = fileparts(outputfilename);
        if length(result)/maxresultsize <= 1
            if isfield(input,'peporglycopep')
                if input.peporglycopep == 1
                    tempoutputfilename = fullfile(p,[f,'_pep','_',num2str(protindse(1)),'_',num2str(protindse(2)),e]);
                    save(fullfile(outputdir,tempoutputfilename),'result','scoreintdata');
                elseif input.peporglycopep == 2
                    tempoutputfilename = fullfile(p,[f,'_N_glycopep','_',num2str(protindse(1)),'_',num2str(protindse(2)),e]);
                    save(fullfile(outputdir,tempoutputfilename),'result','scoreintdata');
                elseif input.peporglycopep == 3
                    tempoutputfilename = fullfile(p,[f,'_O_glycopep','_',num2str(protindse(1)),'_',num2str(protindse(2)),e]);
                    save(fullfile(outputdir,tempoutputfilename),'result','scoreintdata');
                end
            else
                tempoutputfilename = fullfile(p,[f,'_',num2str(protindse(1)),'_',num2str(protindse(2)),e]);
                save(fullfile(outputdir,tempoutputfilename),'result','scoreintdata');
            end
%            save(fullfile(outputdir,tempoutputfilename),'result','scoreintdata');
        else
            resultbackup = result;
            groupsize = ceil(length(resultbackup)/maxresultsize);
            for ii = 1:groupsize
                resultgroupse = [(ii-1)*maxresultsize+1,min(length(resultbackup),ii*maxresultsize)];
                tempoutputfilename = fullfile(p,[f,'_',num2str(protindse(1)),'_',...
                    num2str(protindse(2)),'_Part_',num2str(ii),'of',num2str(groupsize),e]);
                result = resultbackup(resultgroupse(1):resultgroupse(2));
                save(fullfile(outputdir,tempoutputfilename),'result','scoreintdata');
                pause(5);
            end
        end
    case 'browsersave'
        scoreintdata = scoreoptions;
        scoreintdata.filetype = 'BrowserSave';
        scoreoptions = scoreintdata.scoreoptions;
        input = scoreintdata.sliminput;
        [~,f,e] = fileparts(outputfilename);
        tempoutputfilename = [f,e];
    case 'presearchsave'
        prot = allpepdata.glypep;
        protname = allpepdata.FASTAhead;
        sliminput = input;
        fldnms = fieldnames(sliminput);
        if ismember('handles',fldnms)
            sliminput = rmfield(sliminput,'handles');
        end
        if ismember('allfragments',fldnms)
            sliminput = rmfield(sliminput,'allfragments');
        end
        if ismember('allfragments_presearch',fldnms)
            sliminput = rmfield(sliminput,'allfragments_presearch');
        end
        scoreintdata.scoreoptions = scoreoptions;
        scoreintdata.sliminput = sliminput;
        scoreintdata.prot = prot;
        scoreintdata.protname = protname;
        scoreintdata.filetype = 'PreSearchSave';
        if scoreoptions.isHCDtrigger
            SASSO = triggerdata.SASSO;
            colnames = triggerdata.colnames;
            colselected = ismember(upper(colnames),upper(scoreoptions.analyzefragmode));
            colnames = colnames(colselected);
            SASSO = [SASSO(:,colselected),SASSO(:,end-1:end)];
            scoreintdata.SASSO = SASSO;
            scoreintdata.colnames = colnames;
        end
        [p,f,e] = fileparts(outputfilename);
        if length(result)/maxresultsize <= 1
            tempoutputfilename = fullfile(p,[f,'_presearch_',num2str(protindse(1)),'_',num2str(protindse(2)),e]);
            save(fullfile(outputdir,tempoutputfilename),'result','scoreintdata');
        else
            resultbackup = result;
            groupsize = ceil(length(resultbackup)/maxresultsize);
            for ii = 1:ceil(length(resultbackup)/maxresultsize)
                resultgroupse = [(ii-1)*maxresultsize+1,min(length(resultbackup),ii*maxresultsize)];
                tempoutputfilename = fullfile(p,[f,'_presearch_',num2str(protindse(1)),'_',...
                    num2str(protindse(2)),'_Part_',num2str(ii),'of',num2str(groupsize),e]);
                result = resultbackup(resultgroupse(1):resultgroupse(2));
                save(fullfile(outputdir,tempoutputfilename),'result','scoreintdata');
                pause(5);
            end
        end
end
end