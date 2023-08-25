classdef Chemformula < handle
    %CHEMFORMULA: An object representing chemical formula
    %
    %  Syntax;
    %      cfObj = Chemformula(cfstruct);
    %
    %
    %  Example:
    %     H02CFObj  = Chemformula(struct('H',2,'O',1));
    %     HObj      = Chemformula(struct('O',1));
    %     H2O2CFobj = H02CFObj.clone;
    %     H2O2CFobj.add(HObj);
    %
    % See also CHEMELE GLYCAN AMINOACID PROTEASE MODIFICATION
    % Author: Gang Liu
    % Date: 1/19/2015
    % Last updated: Kai Cheng
    % Date: 11/10/2020
    
    properties
        cfstruct
    end
    
    methods 
        % Create an object for inputs
        function obj=Chemformula(varargin)
            if(nargin==1)
                cffieldnames = fieldnames(varargin{1});
                for i = 1:length(cffieldnames)
                    obj.cfstruct.(upper(cffieldnames{i})) = varargin{1}.(cffieldnames{i});
                end
            elseif(nargin==0)
                obj.cfstruct=[];
            else
                error('MATLAB:GlycoPAT:WRONGNUMINPUT','Incorrect number of inputs');
            end
        end
    end
    
    methods
        % Convert whole chemical formula to string, including number of
        %     atoms
        function cfstring = tostring(obj)
            fnames = fieldnames(obj.cfstruct);
            cfstring='';
            for i = 1 : length(fnames)
                cfstring = [cfstring,fnames{i}];
                cfstring = [cfstring,num2str(obj.cfstruct.(fnames{i}))];
            end
        end
        
        % Check if element names are supported
        function isvalidelename = checkelename(obj)
            formulaA = obj.cfstruct;
            eleA = fieldnames(formulaA);
            isvalidelename = true;
            for i = 1:length(eleA)
                if(sum(strcmpi(eleA{i},Chemele.valelename))==0)
                    isvalidelename = false;
                    return
                end
            end
        end
        
        % Merge a new formula with existing one 
        function add(obj,formB)
            formulaA = obj.cfstruct;
            formulaB = formB.cfstruct;
            
            if(isempty(formulaA))
                obj.cfstruct=formulaB;
                return;
            elseif(isempty(formulaB))
                obj.cfstruct=formulaA;
                return
            end
            
            eleA = fieldnames(formulaA);
            eleB = fieldnames(formulaB);
            eleC = eleA;
            
            for i=1:length(eleB)
                if(sum(strcmp(eleB{i},eleA))==0)
                    formulaA.(eleB{i})=0;
                    eleC{length(eleC)+1}=eleB{i};
                end
            end
            
            for i=1:length(eleC)
                if(sum(strcmp(eleC{i},eleB))==0)
                    formulaB.(eleC{i})=0;
                end
            end
            
            for i=1:length(eleC)
                obj.cfstruct.(eleC{i})=formulaA.(eleC{i})+formulaB.(eleC{i});
            end
        end
        
        % Remove a formula from an existing one
        function sub(obj,formB)
            formulaA= obj.cfstruct;
            formulaB= formB.cfstruct;
            
            if(isempty(formulaA))
                eleB = fieldnames(formulaB);
                for i = 1 : length(eleB)
                    obj.cfstruct.(eleB{i})= formulaB.(eleB{i})*-1;
                end
                return
            elseif(isempty(formulaB))
                obj.cfstruct=formulaA;
                return
            end
            
            eleA = fieldnames(formulaA);
            eleB = fieldnames(formulaB);
            
            eleC = eleA;
            
            for i=1:length(eleB)
                if(sum(strcmpi(eleB{i},eleA))==0)
                    formulaA.(eleB{i})=0;
                    eleC{length(eleC)+1}=eleB{i};
                end
            end
            
            for i=1:length(eleC)
                if(sum(strcmpi(eleC{i},eleB))==0)
                    formulaB.(eleC{i})=0;
                end
            end
            
            for i=1:length(eleC)
                objcfstruct.(eleC{i})=formulaA.(eleC{i})- formulaB.(eleC{i});
            end
            
            obj.cfstruct=objcfstruct;
        end
        
        % Multiply the number of atoms in an existing formula
        function multiply(obj,numberStruct)
            if(numberStruct==0)
                error('MATLAB:GlycoPAT:WRONGINPUT','Input must not be equal to 0')
            end
            
            formulaA=obj.cfstruct;
            
            eleA=fieldnames(formulaA);
            for i=1:length(eleA)
                objcfstruct.(eleA{i}) = ...
                    formulaA.(eleA{i}) * ...
                    numberStruct;
            end
            
            obj.cfstruct=objcfstruct;
        end
        
        % Duplicate an existing formula
        function newobj=clone(obj)
            newobj=Chemformula;
            newobj.cfstruct = obj.cfstruct;
        end
    end
end