% Main script for running the surface - subsurface model
% Here you can define which year you want to compute, define your
% parameters, plot output variables, evaluate model performance....
%
% Author: Baptiste Vandecrux (bava@byg.dtu.dk)
% ========================================================================
clearvars
close all
clc
addpath(genpath('.\lib'))
addpath(genpath('Input'),genpath('Output'))

% All parameters are defined in a csv files in Input folder, however it is
% possible to overwrite the default value by defining them again in the
% "param" struct hereunder.

% Increasing resolution towards the sufrace
% param.cdel = zeros(32,1);
% param.cdel = 3.^([1:32]'-29)+0.15;

% High resolution grid (comment if not needed)
NumLayer = 200;
param.z_max = 40;
param.dz_ice = param.z_max/NumLayer;
param.verbose = 1;
param.lim_new_lay = 0.04;

% param.ConductionModel = 0;      % if 1 does CONDUCTION ONLY
% In the conduction model, the temperature profile is reseted every night
% at 2am local time usingg thermistor string readings

% Heterogeneous precolation from Marchenko et al. (2017)
% this considers only redistribution of the water from the first layer into
% the rest of the subsurface
% An alternative is to go through the whole column and check whether piping
% can occur at any depth
param.hetero_percol = 0; % 1 = whole scheme on; 0 = standard percolation
param.hetero_percol_p = 1; % binomial probability for a piping event to be initiated
                            % In Marchenko et al. (2017) this happens at
                            % every time step (probability 1)
param.hetero_percol_frac = 1; % fraction of the available water that can be
                            % In Marchenko et al. (2017) all the available
                            % water goes into redistribution (frac = 1)
param.hetero_percol_dist = 5; % characteristic distance until which
                            % preferential percolation operates
                            % When using uniform probability distribution 
                            % for the redistribution Marchenko et al. (2017)
                            % recommends between 4.5 and 6 m as cut-off value

% result of a tuning of densification schemes
%    param.a_dens = 30.25;
%    param.b_dens = 0.7;

param.year    = 0;
% by defining param.year, the model will be run only for that melt year
% (i.e. 1st. april to 1st april) however you still need to make sure that 
% "rows"  is set so that the appropriate values will be read in the weather 
% data. To run the model only from 1 to rows, just leave param.year=0.

% param.track_density = 0;
param.track_temp = 0;
param.avoid_runoff = 1; % 0 = runoff calculated with Zuo & Oerlemans (1996)

param.retmip = 0;

station_list = {'DYE-2','CP1', 'Summit','NASA-SE',...
     'NASA-E','NASA-U','TUNU-N','SouthDome','Saddle'};
%{'KAN-U','NGRIP','Saddle'};
% {'KAN-U'};
% {'Dye-2_16','Dye-2_long','Summit','FA'};

for i =1:length(station_list)
    param.station =  station_list{i}; 
%     param.InputAWSFile = sprintf(['U:\\Storage Baptiste\\Code\\FirnModel_bv_v1.3\\Input' ...
%         '\\Weather data\\Temperature tracking\\data_%s_combined_hour.txt'],param.station);

switch param.station
    case 'CP1'
%         param.InputAWSFile = 'data_CP1_1998-2010.txt';
        param.InputAWSFile = '../../Data/AWS/Output/CP1/data_CP1_combined_hour.txt';
    case 'DYE-2'
%         param.InputAWSFile = 'data_DYE-2_1998-2015.txt';
%         param.InputAWSFile = 'data_DYE-2_combined_hour.txt';
%         param.InputAWSFile = 'data_DYE-2_combined_hour_cor.txt';
%         param.InputAWSFile = 'data_DYE-2_restricted_hour.txt';
        param.InputAWSFile = '../../Data/AWS/Output/DYE-2/data_DYE-2_combined_hour.txt';
    case {'DYE-2_long' 'Dye-2_long'}
        param.InputAWSFile = 'data_DYE-2_1998-2015.txt';
    case {'DYE-2_HQ' 'Dye-2_16'}
        param.InputAWSFile = 'data_DYE-2_Samira_hour.txt';

    case 'Summit'
%         param.InputAWSFile = 'data_Summit_1990-2015.txt';
%         param.InputAWSFile = 'data_Summit_2000-2015.txt';
%         param.InputAWSFile = 'data_Summit_combined_hour.txt';
        param.InputAWSFile = '../../Data/AWS/Output/Summit/data_Summit_combined_hour.txt';
    case 'NASA-SE'
%         param.InputAWSFile = 'data_NASA-SE_1998-2015.txt';
        param.InputAWSFile = '../../Data/AWS/Output/NASA-SE/data_NASA-SE_combined_hour.txt';
    case 'NASA-E'
        param.InputAWSFile = '../../Data/AWS/Output/NASA-E/data_NASA-E_combined_hour.txt';
    case 'NASA-U'
        param.InputAWSFile = '../../Data/AWS/Output/NASA-U/data_NASA-U_combined_hour.txt';
    case 'SouthDome'
        param.InputAWSFile = '../../Data/AWS/Output/SouthDome/data_SouthDome_combined_hour.txt';
    case 'TUNU-N'
        param.InputAWSFile = '../../Data/AWS/Output/TUNU-N/data_TUNU-N_combined_hour.txt'; 
  case 'Saddle'
        param.InputAWSFile = '../../Data/AWS/Output/Saddle/data_Saddle_combined_hour.txt'; 
  case 'NGRIP'
        param.InputAWSFile = '../../Data/AWS/Output/NGRIP/data_NGRIP_combined_hour.txt'; 
        
% ======== other stations =================
    case {'KAN-U' 'KAN_U'}
%         param.InputAWSFile = 'data_KAN_U_combined_hour.txt';
        param.InputAWSFile = 'data_KAN_U_2012.txt';
%         param.InputAWSFile = '../../Data/AWS/Output/KAN_U/data_KAN_U_combined_hour.txt';
    case {'Miege' 'FA'}
        param.InputAWSFile = 'data_Miege_combined_hour.txt';
    case 'NUK_K'
        param.InputAWSFile = 'data_NUK_K_combined_hour.txt';
        param.shallow_firn_lim = 3;
        param.deep_firn_lim = 5;
        param.min_tot_thick = 15;

        param.lim_new_lay = 0.005;
        param.z_max = 20;
        param.dz_ice = param.z_max/NumLayer;
    otherwise
        disp('Missing data file for the requested station.');
end

% When you add sites
% 1) in Main.m: define path of InputAWSFile
% 2) in  .m: Define path for density profile
% 3) Make sure all information is reported in Input/Constants/InfoStation.csv file

%%  ========= Run model with name tag of your choice =======================
[RunName, c] = HHsubsurf(param);
end
% 