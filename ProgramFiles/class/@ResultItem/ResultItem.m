classdef ResultItem
    % RESULTITEM: An object storing the properties of all fields in result
    %
    % Syntax:
    % allresultfieldnames = ResultItem.allresultitems
    % numericfields = ResultItem.itemisnumeric
    % numcellfields = ResultItem.itemisnumcell
    % stringfields = ResultItem.itemisstring
    %
    % Input:
    % N/A
    %
    % Output;
    % allresultfieldnames, numericfields, numcellfields,
    %     stringfields: 1 x n cell array of strings. The fieldnames that are of
    %     a specific type.
    %
    % Examples:
    %     allresultfieldnames = ResultItem.allresultitems
    %
    % allresultfieldnames =
    %
    %   1×29 cell array
    %
    %   Columns 1 through 7
    %
    %     {'Scan'}    {'Mono'}    {'Theo'}    {'Expt'}    {'Charge'}    {'SGP'}    {'Fragmode'}
    %
    %   Columns 8 through 13
    %
    %     {'PeakLag'}    {'HtCenter'}    {'HtAvg'}    {'PercentIonMatch'}    {'Pvalue'}    {'DecoyRatio'}
    %
    %   Columns 14 through 19
    %
    %     {'Top10'}    {'SelectPeak'}    {'NpFrag'}    {'NgFrag'}    {'NmFrag'}    {'Enscore'}
    %
    %   Columns 20 through 25
    %
    %     {'DecoyEnscore'}    {'Retime'}    {'Quant'}    {'Glycov'}    {'Y0Y1Y2'}    {'ProteinID'}
    %
    %   Columns 26 through 29
    %
    %     {'Fracpepfragmatch'}    {'Fracglyfragmatch'}    {'MDiff'}    {'Protein'}
    %
    % Note:
    % N/A
    %
    % See also: CHEMELE  GLYCAN  CHEMFORMULA  MODIFICATION  AMINOACID
    %     RESULTITEM
    
    % Author: Kai Cheng
    % Date:  1/19/2020
    % Date Lastly Updated: 11/10/2020
    
    properties (Constant)
        % All field names in the "header" data structures, this structure
        %     stores scoring parameters
        allheaderitems = {'PeptideFile','ExperimentData','IsHCDTrigger',...
            'MS1Tolerence','MS2Tolerence','FragmentationNumber',...
            'LabileMonosaccharide','GlycanPeptideSimultaneousFragmentation','AdditionalOxoniumIon',...
            'PeptideFragmentType','MaximumGlycanStubLength','CutOffMedian',...
            'FracMax','SelectivePeak'};
        
        % All field names in the "result" data structures, this structure
        %     stores scoring results
        allresultitems = {'Scan','Mono','Theo',...
            'Expt','Charge','SGP',...
            'Fragmode','PeakLag','HtCenter',...
            'HtAvg','PercentIonMatch','Pvalue',...
            'DecoyRatio','Top10','SelectPeak',...
            'NpFrag','NgFrag','NmFrag',...
            'Enscore','DecoyEnscore','Retime',...
            'Quant','Glycov','Pepcov',...
            'Y0Y1Y2','ProteinID','Fracpepfragmatch',...
            'Fracglyfragmatch','MDiff','Protein'};
        
        % Data in "result" that are pure numbers
        itemisnumeric = {'Scan','Mono','Theo',...
            'Expt','Charge','PeakLag',...
            'HtCenter','HtAvg','PercentIonMatch',...
            'Pvalue','DecoyRatio','Top10',...
            'NpFrag','NgFrag','NmFrag',...
            'Enscore','DecoyEnscore','Retime',...
            'Quant','Glycov','Pepcov',...
            'Fracpepfragmatch','Fracglyfragmatch','MDiff'};
        
        % Data in "result" that are cell array of strings converted from
        %     numbers. These data should be converted to numbers before
        %     use
        itemisnumcell = {'SelectPeak','Y0Y1Y2','ProteinID'};
        
        % Data in "result" that are cell array of strings
        itemisstring = {'SGP','Fragmode','Protein'};
    end
end