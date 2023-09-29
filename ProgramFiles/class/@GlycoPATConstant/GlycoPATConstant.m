classdef GlycoPATConstant
    % GLYCOPATCONSTANT: An object storing the information used throughout
    %     GlycoPAT
    %
    % Syntax:
    % combiweight_CID = GlycoPATConstant.combiweight.CID;
    % combiweight_HCD = GlycoPATConstant.combiweight.HCD;
    % combiweight_ETHCD = GlycoPATConstant.combiweight.ETHCD;
    %
    % Input:
    % N/A
    %
    % Output;
    % combiweight: the weight for calculating combiES.
    %
    % Examples:
    % N/A
    %
    % Note:
    % N/A
    
    %  Author: Kai Cheng
    %  Date Lastly Updated: 1/19/2020
    
    properties (Constant)
        % Weight ratio for calculating combined ensemble score (CombiES)
        %     from ES's of different fragmentation methods
        
        combiweight = struct('CID',0.25,...
            'HCD',0.35,...
            'ETHCD',0.4);
        
        % Ensemble score cut-off values for:
        EScutoff_analyzepresearchdata = 0.4;  % In PreSearch, min ES required for further calculation
        EScutoff_fbi = 0.4;  % In FindBestIsomer, min ES required for further calculation
        %Lowered to 0.4 from 0.5 on 12/26/2022
        %03/05/2023 raised from 0.2 to 0.4
        
        % Max. lag of cross-correlation
        Maxlag = 25;
        
        % number of decoy for ES score calulation
        ndecoy = 10;
        
        % Peptide coverage defined as 
        Pepcovcutoff_analyzepresearchdata = 0.3;
        Pepcovcutoff_fbi = 0.0;   % should change to 0.2 afterwards
        Pepcovcutoff_Oglysearch = 0.2; %Lowered to 0.2 on 12/08/2022
        
        TopNProt_Method5_analyzepresearchdata = 2;
        
        KeepTopOglynum = 10000;
        
        fbi_nummmindiagionfragmode = 1;
        fbi_diagnionlogic = 1;%1=and, 2=or
        OGlyPepFDR = 0.01;
        
        N_Proceed_Override = true;
        O_Proceed_Override = true;
        
        
        Scoreweight_HCD_Top10 = 1/6;
        Scoreweight_HCD_Pscore = 0;
        Scoreweight_HCD_Y012 = 1/3;
        Scoreweight_HCD_Glycov = 0;
        Scoreweight_HCD_Pepcov = 1/3;
        Scoreweight_HCD_XCorr = 0;
        Scoreweight_HCD_Oxo = 1/6;
        Scoreweight_HCD_PerIonMatch = 0;
        
        Scoreweight_CID_Top10 = 1/4;
        Scoreweight_CID_Pscore = 1/4;
        Scoreweight_CID_Glycov = 1/4;
        Scoreweight_CID_XCorr = 1/4;
        Scoreweight_CID_Y012 = 0;
        Scoreweight_CID_Pepcov = 0;
        Scoreweight_CID_Oxo = 0;
        Scoreweight_CID_PerIonMatch = 0;
        
        Scoreweight_ETD_Top10 = 0;
        Scoreweight_ETD_Pscore = 1/3;
        Scoreweight_ETD_Glycov = 0;
        Scoreweight_ETD_XCorr = 1/3;
        Scoreweight_ETD_Y012 = 0;
        Scoreweight_ETD_Pepcov = 1/3;
        Scoreweight_ETD_Oxo = 0;
        Scoreweight_ETD_PerIonMatch = 0;
        
        Scoreweight_EThcD_Top10 = 0.25;
        Scoreweight_EThcD_Pscore = 0;
        Scoreweight_EThcD_Glycov = 0.1;
        Scoreweight_EThcD_XCorr = 0.1;
        Scoreweight_EThcD_Y012 = 0;
        Scoreweight_EThcD_Pepcov = 0.45;
        Scoreweight_EThcD_Oxo = 0.1;
        Scoreweight_EThcD_PerIonMatch = 0;
        
        Scoreweight_ETciD_Top10 = 0.1;
        Scoreweight_ETciD_Pscore = 0.7;
        Scoreweight_ETciD_Glycov = 0;
        Scoreweight_ETciD_XCorr = 0.2;
        Scoreweight_ETciD_Y012 = 0;
        Scoreweight_ETciD_Pepcov = 0;
        Scoreweight_ETciD_Oxo = 0;
        Scoreweight_ETciD_PerIonMatch = 0;
        
    end
end