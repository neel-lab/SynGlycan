function output = drawglycan(allsgp,varargin)
% DRAWGLYCAN: draw glycan or glycopeptide. This is the main program
%
% Syntax:
% output = drawglycan(allsgp)
% output = drawglycan(allsgp,options)
% output = drawglycan(allsgp,optionname1,optionvalue1,...)
%
% Input:
% allsgp: string. Glycan and glycopeptide string in IUPAC format. Local
%     options are integrated into this sequence.
% options: structure. User defined global options listed below in Table format.
% optionname and optionvalue. User defined option name and value pairs.
%
% Output:
% output: structure. As of now, it has a field "tcp", an 2 x 2 numerical
% array. First row is the coordinate of figure's lower left corner, second
% row is its upper right corner.
%
% Note: Available global options are (case sensitive):
% ---------------------------------------------------------------------------------------------------------------------------------------------------
%     (Below when the word "normalization" is used, the normalization factor is the distance between the center of
%         two adjacent (connected by a glycosidic bond) monosaccharides. Symbols of monosaccharides will be referred
%         to as "symbols". Lines describing glycosidic bonds will be referred to as "bonds".)
%     Variable name                  Format                                 Description
%     orientation                        String                                    Orientation of glycan structure to be drawn. Available options
%                                                                                                     are: "up", "down", "left" and "right". 
%                                                                                                     Default is "left".  
%     monosacsize                    Double                                  Normalized diameter of symbols.
%                                                                                                     Default is 0.5.
%     msperiwidth                     Double                                 Width of the black rim of symbols in points, 1 point = 1/72 of an
%                                                                                                      inch.
%                                                                                                     Default is 1.
%     bondwidth                        Double                                  Width of bonds in points.
%                                                                                                     Default is 2.
%     perpendicular                  1 x n cell array of                Monosaccharides that will be placed perpendicularly from the
%         monosac                           strings                                   general orientation, the clockwise one will be used first. 
%                                                                                                     Additional symbols will be placed on both sides and distributed
%                                                                                                     evenly.
%                                                                                                     Default is: "Fuc", "Xyl", "Deoxyhexose" and "Pentose".
%     linkinfotheta                     Double                                    C       D        A and B represent symbols. C and D represent the 2 
%     linkinfodist                                                                         /           \           parts of linkage information.
%                                                                                                 A---------B
%                                                                                                 "linkinfotheta" is the angle between AC and AB and the angle
%                                                                                                     between BA and BD. "linkinfodist" is the normalized length of
%                                                                                                     line AC and BD. The linkage information is placed either to the
%                                                                                                     left (when "orientation" is "up" or "down") or below (when 
%                                                                                                     "orientation" is "left" or "right") the bonds. 
%                                                                                                     Default for "linkinfotheta" is -30 (degrees).
%                                                                                                     Default for "linkinfodist" is 0.7.
%     bondbreaksiglength       Double                                 Normalized length of bond cleavage markers.
%                                                                                                     Default is 0.5 (same as "monosacsize").
%     showlink                            Logical                                  Show linkage information or not. 
%                                                                                                     Default is true.
%     fontsize                              Double                                 "General" size of texts (except linkage information). The actual 
%                                                                                                     size of texts is calculated by multiplying a "fontsize" and a
%                                                                                                     factor. These factors are: 
%                                                                                                     1.5: for amino acids,
%                                                                                                             monosaccharides displayed in text form by using local
%                                                                                                                 option "-CHAR",
%                                                                                                             non-glycan PTMs on peptide,
%                                                                                                             first letter of the name of monosaccharide currently does 
%                                                                                                                 not have symbol (represented by white hexagon with
%                                                                                                                 first letter at its center),
%                                                                                                             text realted to anomers displayed by using local options
%                                                                                                                 "-ANOMER";
%                                                                                                     1.25: for monosaccharide information displayed by using local 
%                                                                                                                     options "-U", "-D" or "C",
%                                                                                                                 text related to starting and ending of repeated subunit
%                                                                                                                     displayed by "-RS" and "-RE";
%                                                                                                     1.0: for fragmentation text displayed by "-R", "-NR", "N" and
%                                                                                                                     "C",
%                                                                                                              text related to adduct displayed by "-ADDUCT".
%                                                                                                     Default is 12.
%     workingmode                  String                                   Type of object to be drawn, can be "G" (glycan), "GP" 
%                                                                                                     (glycopeptide) or "P" (peptide). Program will try to decide this
%                                                                                                     by examining input.
%                                                                                                     Default is empty.
%     structspacing                   Double                                 Normalized distance between different glycan structures when
%                                                                                                     there are more than one present.
%                                                                                                     Default is 1.
%     aaspacing                          Double                                 Normalized distance between amino acid letters.
%                                                                                                     Default is 0.75.
%     fileout                                 String                                    The full name of the file where the program will save the drawn
%                                                                                                     figure drawn to. If is empty, the figure will be drawn, but user
%                                                                                                     must manually save it to a file.
%                                                                                                     Default is empty.
%     visible                                 String                                    Controls whether the figure drawn will be displayed. Available
%                                                                                                     options are "yes" and "no".
%                                                                                                     Default is "yes".
%     inputformat                      String                                    The format of input sequence. Currently the following formats
%                                                                                                     are supported: IUPAC, SGP1 and SGP2.
%                                                                                                     Default is "IUPAC"
%     linkfontsize                       Double                                  Font size of linkage information.
%                                                                                                     Default is 12.
%     sortbranches                   Logical                                   If true, when there are more than one monosaccharides attach
%                                                                                                     to a monosaccharide ("branches"), the branches will be
%                                                                                                     sorted based on their linkage information (anomeric
%                                                                                                     carbon (alpha, beta,...) then hydroxyl groups (1,2,...)).
%                                                                                                     Default is true.
%     figurehandle                     Axes object                         The axes where the figure will be drawn upon. If is empty,
%                                                                                                     program will create a figure object and an axes object.
%                                                                                                     Default is empty.
%     resize                                  Logical                                  If true, the size of the figure drawn will be readjusted according
%                                                                                                     to the size of the structure.
%                                                                                                     Default is true;
%     pointsperunit                   Double                                  An constant for estimating points (as defined in "msperiwidth")
%                                                                                                     per normalized length unit.
%                                                                                                     Default is 50.
%     specialoptions                 Structure                               This structure can be used to provide local options to program
%                                                                                                     without encoding them into sequence. To do so, use the
%                                                                                                     name of local options as fieldnames and formatted cell array as
%                                                                                                     values. All local options can be integrated here. The format is:
%                                                                                                     n x 2 cell array, 1st column contains the option values, 2nd 
%                                                                                                     column contains the serial number of the monosaccharide
%                                                                                                     where the options will be applied upon. If an option is applied
%                                                                                                     to more than 1 monosaccharide, place them in different cells.
%                                                                                                     Note: when presenting peptide fragmentations, use "PEPC" and
%                                                                                                     "PEPN" instead of "C" and "N". This is to avoid confusion with 
%                                                                                                     "C" which modifies symbols.
%                                                                                                     Default is empty.
% ---------------------------------------------------------------------------------------------------------------------------------------------------
% For local options available, see GETSEQOPTIONS.
% 
% Example:
% drawglycan('A[GalNAc(a1-3)](-N "B2" -C  "X4")AB[Xyl(b1-3)GalNAc(a1-3)]B(-C "X2.5")C[Fuc(a1-3)]CDD')
%
% Children function:
% USG, CALCGLYPOS, ESTIFIGSIZE, PAINTGLYCAN, GETFRAGINFO, REARRANGEGLYPEP,
% PLOTPEPBREAK
%

% v2.0 update note: additional local options added, see "DrawGlycanPara"
% for details.
% v1.1 update note: draw one structure at a time, the function of drawing
% multiple structure from a cell input has been cancelled.
%

% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2020, Research Foundation for State University of New York. All rights reserved
%

%% INITIATING - BUILD DEFAULT OPTIONS
defoptname = {'orientation','monosacsize','msperiwidth',...
    'bondwidth','perpendicularmonosac',...
    'linkinfodist','linkinfotheta','bondbreaksiglength',...
    'showlink','fontsize','workingmode',...
    'structspacing','aaspacing','fileout',...
    'visible','inputformat','specialoptions',...
    'linkfontsize','sortbranches','figurehandle',...
    'resize','pointsperunit','centertextsize',...
    'centerfontweight','originpointoffset',...
    'zoomfactor'};  % option names

defoptvalue = {'left',.5,1,...
    2,{'Fuc','Xyl','Deoxyhexose','Pentose'},...
    .7,-30,.5,...
    true,12,'',...
    1,.75,'',...
    'on','iupac',{},...
    12,true,[],...
    true,50,11,...
    'normal',[0,0],...
    1};  % option default values

defopt = usg('gen',defoptname,defoptvalue);  % defopt is the option structure, all values are now default
if ~isempty(varargin)
    if isstruct(varargin{1})
        % update option values (input is a structure)
        options = usg('mod',defopt,varargin{1});
    elseif ischar(varargin{1})
        % update option values (input is "option name" - "option value" pairs)
        options = usg('mod',defopt,varargin(1:2:end),varargin(2:2:end));
    elseif iscell(varargin{1})
        % update option values (input is a cell array: {option1,value1,option2,value2,...})
        inputoption = varargin{1};
        options = usg('mod',defopt,inputoption(1:2:end),inputoption(2:2:end));
    end
else
    options = defopt;
end

if ismember(lower(options.orientation),{'right','down'})
    % program is designed for "up" and "left", here the angles need to
    % be processed to avoid mirroring
    options.linkinfotheta = -options.linkinfotheta;
end
options.isaxes = true;  % assuming we have the whole window for us to draw on
% otherwise it's a small area in the window

% glyseq is cell, pepseq is char
if isempty(options.workingmode)
    % split glycan and peptide from input sequence
    [glyseq,pepseq,glypos,~] = distggp(allsgp,options.inputformat);
    if isempty(pepseq)
        options.workingmode = 'G';
    else
        if ~isempty(glyseq)
            options.workingmode = 'GP';
        else
            options.workingmode = 'P';
        end
    end
end
%% handle different input format, produce a "glysgp" and "specialoptions"
glyoptions = [DrawGlycanPara.intglymodinfo,DrawGlycanPara.intglybondinfo,...
    DrawGlycanPara.extglymodinfo,DrawGlycanPara.glyidentityinfo,'CURLYBRACKET'];

if strcmpi(options.inputformat,'IUPAC')
    switch options.workingmode
        case 'G'
            [cleanseq,specialoptions] = getglydressing(strtrim(glyseq));
            glysgp = IUPAC2Sgp(cleanseq,2);
        case 'GP'
            glyseq = cellfun(@strtrim,glyseq,'UniformOutput',false);
            [cleanseq,specialoptions] = cellfun(@getglydressing,glyseq,'uniformoutput',false);
            glysgp = cellfun(@IUPAC2Sgp,cleanseq,num2cell(repmat(2,numel(glyseq),1)),'uniformoutput',false);
            if ~isempty(options.specialoptions)
                extspecialoptions = options.specialoptions;
                extspfldnms = fieldnames(extspecialoptions);
                intspfldnms = fieldnames(specialoptions);
                spfieldnames = setdiff(extspfldnms,intspfldnms);
                for ii = 1:length(spfieldnames)
                    if strcmpi(spfieldnames{ii},'PEPC') || strcmpi(spfieldnames{ii},'PEPN')
                        specialoptions.(spfieldnames{ii}(4)) = extspecialoptions.(spfieldnames{ii});
                    else
                        specialoptions.(spfieldnames{ii}) = [specialoptions.(spfieldnames{ii});...
                            extspecialoptions.(spfieldnames{ii})];
                    end
                end
            end
        case 'P'
            specialoptions = options.specialoptions;
    end
    if ~isempty(specialoptions)  % Program is designed based on SGP format. The monosac. 
        % are arranged differently. Therefore, option's position needs to be adjusted
        if iscell(specialoptions)
            % This means there are multiple glycans to be drawn
            for ii = 1:length(specialoptions)
                if ~isempty(specialoptions{ii})
                    spoptfld = fieldnames(specialoptions{ii});  % special option fields
                    numMono = length(strfind(glysgp{ii},'{'));
                    for jj = 1:length(spoptfld)
                        if ismember(upper(spoptfld{jj}),glyoptions)
                            tempspopt = specialoptions{ii}.(spoptfld{jj});
                            for kk = 1:size(tempspopt,1)
                                tempspopt{kk,2} = bsxfun(@minus,numMono + 1,tempspopt{kk,2});  % convert option associated monosac location from IUPAC to SGP
                            end
                            specialoptions{ii}.(spoptfld{jj}) = tempspopt;
                        end
                    end
                end
            end
        elseif isstruct(specialoptions)
            % Only 1 glycan
            spoptfld = fieldnames(specialoptions);  % special option fields
            numMono = length(strfind(glysgp,'{'));
            for jj = 1:length(spoptfld)
                if ismember(upper(spoptfld{jj}),glyoptions)
                    tempspopt = specialoptions.(spoptfld{jj});
                    for kk = 1:size(tempspopt,1)
                        tempspopt{kk,2} = bsxfun(@minus,numMono + 1,tempspopt{kk,2});  % convert option associated monosac location from IUPAC to SGP
                    end
                    specialoptions.(spoptfld{jj}) = tempspopt;
                end
            end
        end
    else
        specialoptions = options.specialoptions;
    end
elseif strcmpi(options.inputformat,'GLYCAM')
    % NOT FINISHED YET - STILL IN TEST PHASE
    if ~any(strfind(allsgp,'('))
        allsgp = glycam2iupac(allsgp);  % not tested
    end
    options.inputformat = 'IUPAC';
    [cleanseq,specialoptions] = getglydressing(strtrim(allsgp));
    glysgp = IUPAC2Sgp(cleanseq,2);
elseif strcmpi(options.inputformat,'SGP1')
    % Specially designed for GlycoPAT 2 use - "specialoptions" will be used
    %     instead of local options.
    % Assuing user has already defined "options.workingmode".
    switch upper(options.workingmode)
        case 'G'
%             glysgp = sgp1to2(allsgp,1);
            glysgp = sgp1to2a(allsgp,1);
            specialoptions = options.specialoptions;
        case 'GP'
            [glyseq,pepseq,glypos,PTMid] = distggp(allsgp,options.inputformat);
            % "glyseq" contains non-glycan PTMs, for consistency in
            % programming, put non-glycan PTMs back to peptide sequence.
            % They are later handled by a field "NGMOD"
            glyseq = cellfun(@strtrim,glyseq,'UniformOutput',false);
            glysgp = cellfun(@sgp1to2a,glyseq,num2cell(PTMid),'UniformOutput',false); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%removed second input 
            %             glysgp = cellfun(@sgp1to2b,glyseq,num2cell(PTMid),'UniformOutput',false);
            specialoptions = options.specialoptions;
            if isempty(specialoptions)
                specialoptions = cell(size(glysgp));
            end
        case 'P'
            pepseq = allsgp;
            specialoptions = options.specialoptions;
    end
elseif strcmpi(options.inputformat,'SGP2')
    glysgp = allsgp;
    specialoptions = options.specialoptions;
elseif strcmpi(options.inputformat,'LINUCS')
    % TO BE COMPLETED
elseif strcmpi(options.inputformat,'BCSDB')
    % TO BE COMPLETED
elseif strcmpi(options.inputformat,'LINEARCODE')
    % TO BE COMPLETED
elseif strcmpi(options.inputformat,'GLYCOCT')
    % TO BE COMPLETED
    % sequence itself is a connection table, reformat and put it in specialoptions.
elseif strcmpi(options.inputformat,'WURCS')
    % TO BE COMPLETED
elseif strcmpi(options.inputformat,'KCF')
    % TO BE COMPLETED
elseif strcmpi(options.inputformat,'CABOSML')
    % TO BE COMPLETED
elseif strcmpi(options.inputformat,'GLYDE')
    % TO BE COMPLETED
else
    return
end
if isempty(options.figurehandle)
    parentfigure = figure('visible',options.visible);
    options.figurehandle = axes(parentfigure,'Units','normalized','Position',[0 0 1 1],'Color','None');
    hold on;
    options.isaxes = false;
end
ax = options.figurehandle;
switch upper(options.workingmode)
    case 'G'
        [plotinfoout,specialoptions] = calcglypos(glysgp,options,specialoptions);  % only draw 1 glycan structure, for multiple structures, call func. multiple times
        plotinfoout.mspos = plotinfoout.mspos + options.originpointoffset;
        % specialoptions are merged into plotinfoout
        if ~isempty(specialoptions)
            plotinfoout = calcattachmentpos(plotinfoout,options);
        end
        [plotinfoout,glytcp] = paintglycan(plotinfoout,options,[]);
        atmtcp = paintattachment(plotinfoout,options,glytcp);
        output.tcp = [min(glytcp(1,1),atmtcp(1,1)),min(glytcp(1,2),atmtcp(1,2));...
            max(glytcp(2,1),atmtcp(2,1)),max(glytcp(2,2),atmtcp(2,2))];
        estimatefigsize(plotinfoout,'',output.tcp,options);
    case 'GP'  % default format: IUPAC
        options.orientation = 'up';
        plotinfoout = cell(size(glysgp));
        for ii = 1:length(glysgp)
            [plotinfoout{ii},specialoptions{ii}] = calcglypos(glysgp{ii},options,specialoptions{ii});  % only draw 1 glycan structure, for multiple structures, call func. multiple times
            % specialoptions are merged into plotinfoout
            if ~isempty(specialoptions{ii})
                plotinfoout{ii} = calcattachmentpos(plotinfoout{ii},options);
            end
        end
        [pepoptions,cleanpepseq] = getseqoptions(pepseq,'p');
        if ~isempty(specialoptions{1}) && ismember('pepC',fieldnames(specialoptions{1}))
            if isempty(pepoptions) || ~isfield(pepoptions,'C')
                pepoptions.C = specialoptions{1}.pepC;
            else
                pepoptions.C = [pepoptions.C;specialoptions{1}.pepC];
            end
        end
        if ~isempty(specialoptions{1}) && ismember('pepN',fieldnames(specialoptions{1}))
            if isempty(pepoptions) || ~isfield(pepoptions,'N')
                pepoptions.N = specialoptions{1}.pepN;
            else
                pepoptions.N = [pepoptions.N;specialoptions{1}.pepN];
            end
        end
        
        [plotinfoout,peppos] = rearrangeglypep(cleanpepseq,glypos,plotinfoout,options);
        pepobj = cell(1,length(cleanpepseq));
        for ii = 1:length(cleanpepseq)
            pepobj{ii} = text(ax,peppos(ii),-1,char(double(cleanpepseq(ii))),'FontName','Helvetica','VerticalAlignment','top',...
                'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
        end
        tcp = zeros(2,2);
        tcp(1,1) = pepobj{1}.Extent(1);
        tcp(1,2) = pepobj{1}.Extent(2);
        tcp(2,1) = peppos(end) + options.aaspacing;
        glytcp = tcp;
        atmtcp = tcp;
        for ii = 1:length(plotinfoout)
            if ~isempty(plotinfoout{ii})
                [plotinfoout{ii},glytcp] = paintglycan(plotinfoout{ii},options,glytcp);
                atmtcp = paintattachment(plotinfoout{ii},options,atmtcp);
            end
        end
    case 'P'
        [pepoptions,cleanpepseq] = getseqoptions(pepseq,'p');
        peppos = options.aaspacing * [0:length(cleanpepseq)-1];
        tcp = zeros(2);
        tcp(1,1) = peppos(1) - options.aaspacing;
        tcp(1,2) = -2;
        tcp(2,1) = peppos(end) + options.aaspacing;
        glytcp = zeros(2);
        atmtcp = zeros(2);
        if ~isempty(specialoptions) && ismember('pepC',fieldnames(specialoptions{1}))
            pepoptions.C = specialoptions{1}.pepC;
        end
        if  ~isempty(specialoptions) && ismember('pepN',fieldnames(specialoptions{1}))
            pepoptions.N = specialoptions{1}.pepN;
        end
        pepobj = cell(1,length(cleanpepseq));
        for ii = 1:length(cleanpepseq)
            pepobj{ii} = text(ax,peppos(ii),-1,char(double(cleanpepseq(ii))),'FontName','Helvetica','VerticalAlignment','top',...
                'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
        end
        glytcp = zeros(2);
        atmtcp = zeros(2);
        tcp = zeros(2);
end
if ~options.isaxes
    set(ax,'visible','off')
end
axis(options.figurehandle,'equal')
if strcmpi(options.workingmode,'gp') || strcmpi(options.workingmode,'p')
    pbtcp = plotpepbreak(peppos,pepoptions,glypos,options);  % this is the last step because frag marker length is related to figure size
    if isfield(pepoptions,'NGMOD')
        ngmods = pepoptions.NGMOD;
        for ii = 1:size(ngmods,1)
            text(ax,peppos(ngmods{ii,2}),-2,ngmods{ii,1},'FontName','Helvetica','VerticalAlignment','bottom',...
                'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
        end
    end
    tcp = [min([glytcp(1,1),atmtcp(1,1),tcp(1,1),pbtcp(1,1)]),min([glytcp(1,2),atmtcp(1,2),tcp(1,2),pbtcp(1,2)]);...
        max([glytcp(2,1),atmtcp(2,1),tcp(2,1),pbtcp(2,1)]),max([glytcp(2,2),atmtcp(2,2),tcp(2,2),pbtcp(2,2)])];
    output.tcp = tcp;
    if strcmpi(options.workingmode,'gp')
        set(ax,'XLim',[tcp(1,1),tcp(2,1)],'YLim',[tcp(1,2),tcp(2,2)+options.monosacsize]);
    elseif strcmpi(options.workingmode,'p')
        set(ax,'XLim',[peppos(1) - options.aaspacing,peppos(ii) + options.aaspacing],'YLim',[-2.5,1]);
    end
end

if ~isempty(options.fileout)
    set(gcf, 'InvertHardcopy', 'off')
    saveas(gcf,options.fileout)
end
if strcmpi(options.visible,'off')
    p=get(options.figurehandle,'Parent');
    close(p);
end
end