function iso_dist = customIsoDist(chnos,mass_max,mass_int,resDa)
%CUSTOMISODIST: This function performs quantitative analysis on a
% given msdata structure. 
%
% Syntax: 
%   iso_dist = customIsoDist(chnos,mass_max,mass_int,resDa)
%
% Input: 
%   CHNOS vector for composition ratios, maximum mass of composition
%   in distribution, interval between masses, resolution in Da
%
% Output: 
%   isotopic distribution cell array containing distribution with
%   corresponding compositions
%
% Examples:
%   iso_dist = customIsoDist([],10000,1,0.0625)
%
%See also: PREPROCESSGUI, DEFAULTMONOISO, GETDEFAULTMONOISOWDIST,
% QUANTITATIVEANALYSIS
%

f = waitbar(0,'','Name','Calculating Isotopic Distribution...');
% [0.024321959339170,0.476559048284606,0.001461698046926,0.012941719677186,0] GLYCAN CHNOS VALUES
% [0.0411425801242450,0.0667414095012345,0.00760993965077790,0.0203467465280654,0.000219515203030292] GLYCOPEPTIDE CHNOS VALUES
if ~isempty(chnos)
    avggp = chnos;
else
    avggp = [0.0411425801242450,0.0667414095012345,0.00760993965077790,0.0203467465280654,0.000219515203030292]; % [C H N O S]
end
if mass_max == 0
    mass_max = Inf;
end
int_length = ceil(mass_max)/mass_int;
iso_dist = cell(int_length,2);
for ii = 1:int_length
    mass_chnos = mass_int*ii*avggp;
    chem_formula = strcat('C',num2str(round(ii*avggp(1),4)),'H',num2str(round(ii*avggp(2),4)),'N',num2str(round(ii*avggp(3),4)),'O',num2str(round(ii*avggp(4),4)),'S',num2str(round(ii*avggp(5),4)));
    iso_dist{ii,1} = isotopicdist(mass_chnos,'Resolution',resDa,'showplot',false);
    iso_dist{ii,2} = chem_formula;
    waitbar(ii/ceil(mass_max),f,sprintf('%d/%d',ii,ceil(mass_max)));
end
close(f);
end