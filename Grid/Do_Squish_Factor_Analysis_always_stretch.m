function [Lambda_Results] = Do_Squish_Factor_Analysis_always_stretch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will run the expansion_factor script on a spreadsheet of data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spreadsheet{1} = uigetfile('*.xlsx','Choose wt_c Spreadsheet');
spreadsheet{2} = uigetfile('*.xlsx','Choose trip Spreadsheet');
spreadsheet{3} = uigetfile('*.xlsx','Choose wt_e Spreadsheet');

cd 'W:\RMunn\Munn_Mallory_Hardcastle_Chetkovich_Giocomo_Cell_Submission\R_analysis\Grid\Definitive_Lambda\Always_stretch_small_map'

for h = 1:3
% get data 
[~,~,data] = xlsread(spreadsheet{h});
% load the session, tetrode, and unit to load the right data
row1 = data(1,:);
kind = h; %% 1 FOR WT_C, 2 FOR TRIP, 3 FOR WT_E
Session_square = data(2:end,(strcmp(row1, 'Sessions_square') == 1));
Session_rectangle = data(2:end,strcmp(row1, 'Sessions_rectangle') == 1);
Tetrode_square = data(2:end,find(strcmp(row1, 'Tetrode_square') == 1,1,'first'));
Tetrode_rectangle = data(2:end,find(strcmp(row1, 'Tetrode_rectangle') == 1,1,'first'));
Unit_square = data(2:end,(strcmp(row1, 'Unit_square') == 1));
Unit_rectangle = data(2:end,(strcmp(row1, 'Unit_rectangle') == 1));
Lambda_Results = table;

for n = 1:length(Session_square)
    file2 = strcat(Session_square{n},'_pos.mat');
    file1 = strcat(Session_square{n},'_T',num2str(Tetrode_square{n}),'C',num2str(Unit_square{n}),'.mat');
    file4 = strcat(Session_rectangle{n},'_pos.mat');
    file3 = strcat(Session_rectangle{n},'_T',num2str(Tetrode_rectangle{n}),'C',num2str(Unit_rectangle{n}),'.mat');
 
    [rho_max(n,1),mean_rho(n,1),yShift(n,1),xShift(n,1),lambda(n,1)] = squish_factor_translation_always_stretch(file1, file2, file3, file4); 
    counter = sprintf('%s%s%s%s', 'For Cell ',num2str(n),' of ',num2str(length(Session_square)));
    disp(counter);
end

Lambda_Results.rho_max = rho_max; Lambda_Results.rho_mean = mean_rho; Lambda_Results.yShift_cm = yShift; Lambda_Results.xShift_cm = xShift;
Lambda_Results.Lambda = lambda;

if kind == 1
    Lambda_Results_wt_c_stretch = Lambda_Results;
    save wt_c_lambda_results_stretch_rectangle_26cm.mat Lambda_Results_wt_c_stretch
elseif kind == 2
    Lambda_Results_t8b_stretch = Lambda_Results;
    save t8b_lambda_results_stretch_rectangle_26cm.mat Lambda_Results_t8b_stretch
elseif kind == 3
    Lambda_Results_wt_e_stretch = Lambda_Results;
    save wt_e_lambda_results_stretch_rectangle_26cm.mat Lambda_Results_wt_e_stretch
end
clear Lambda_Results
end

