% This script helps to generate a density profile to be used as initial
% state for the subsurface model
%
% Author: Baptiste Vandecrux (bava@byg.dtu.dk)
% ========================================================================

% clear all
% close all
load Core_all

station = 'NGRIP';
% Plot in the site you want
%  PlotCore(Core,'Site',station)

% Find the index of the core you want to use
% CoreList(CoreAvg)
% search by name rather than by index since index can change in the future

switch station
    case 'CP1'
        ind = FindCore(Core,'Name','CORE 6945');
    case 'DYE-2'
        ind = FindCore(Core,'Name','DYE2 1998 core B');
    case 'DYE-2_HQ'
        ind = FindCore(Core,'Name','core_10_2016');
    case 'KAN-U'
        ind = FindCore(Core,'Name','core_1_2012');
    case 'NASA-SE'
        ind = FindCore(Core,'Name','CORE 6642 (B)');
    case 'Summit'
        ind = FindCore(Core,'Name','T99_1990');
%         ind = FindCore(Core,'Name','Dome GRIP');
    case 'Miege'
        ind = 120; %or 8
    case 'NASA-U'
        ind = FindCore(Core,'Name','CORE 7347');
    case 'SouthDome'
        ind = FindCore(Core,'Name','S. Dome Core A');
%         ind = FindCore(Core,'Name','S. Dome Core B');
    case 'NASA-E'
        ind = FindCore(Core,'Name','NASA East Core A');
    case 'Saddle'
        ind = FindCore(Core,'Name','N. Dye 3 (Saddle) - A');
    case 'TUNU-N'
        ind = FindCore(Core,'Name','B19_NGT19_1994');
    case 'NGRIP'
        ind = FindCore(Core,'Name','NG97S2~1-3bag');
%         ind = FindCore(Core,'Name','NGRIP2001S5');
       
end

% Creating density profile and writing it into the Input folder
depth = Core{ind}.Data.Depth/100;
density = Core{ind}.Data.Density;
ice_perc = Core{ind}.Data.Type_perc;

% if strcmp(station,'KAN-U')
%     depth = depth';
%     density = density';
%     ice_perc = ice_perc';
%     ice_perc(length(ice_perc):length(density)) = 0;
% end
% DensProfile = [depth, density, ice_perc];
filename = sprintf('./Input/Initial state/DensityProfile_%s_%i.csv',station,Core{ind}.Info.DateCored.Year);
    M  = [depth, density];
    M_table = array2table(M,'VariableName', {'depth_m', 'density_kgm3'});

    writetable(M_table,filename,'Delimiter',';')

fprintf('Initial density profile was generated from core %s and placed in Input folder.\n',Core{8}.Info.Name)
PlotCore(Core,'CoreNumber',ind);

