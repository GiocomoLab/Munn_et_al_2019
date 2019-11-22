
load('Z:\Users\RMunn\R_analysis\Head Direction\Compression\Model_ID_HD_Cells_Basic_Stats_and_Raw_Tuning_Curves_WT_C.mat');
wt_c_of = Results_hdof;
wt_c_sq = Results_hdsq;

clear Results_hdof Results_hdsq;

load('Z:\Users\RMunn\R_analysis\Head Direction\Expansion\Model_ID_HD_Cells_Basic_Stats_and_Raw_Tuning_Curves_WT_E.mat');
wt_e_of = Results_hdof;
wt_e_sq = Results_hdsq;

wt_c_mean_fr_diff = wt_c_sq.meanFR_sq - wt_c_of.meanFR;
wt_e_mean_fr_diff = wt_e_sq.meanFR_sq - wt_e_of.meanFR;


for j = 1:height(wt_c_of)
    wt_c_tuning_curve_of(j,:) = wt_c_of.Raw_Curve{j,1};
    wt_c_tuning_curve_sq(j,:) = wt_c_sq.Raw_Curve{j,1};
end
hd_tuning_axis = linspace(0,2*pi,60);

wt_c_sq_max = (max(wt_c_tuning_curve_sq,[],2));
wt_c_of_max = (max(wt_c_tuning_curve_of,[],2));

for l = 1:size(wt_c_sq_max,1)
[h,hg,ia] = intersect(wt_c_sq_max(l),wt_c_tuning_curve_sq(l,:),'stable');
idxwt(l,1) = ia;
pref_ang_sq(l,1) = hd_tuning_axis(idxwt(l));
[h,hg,ia] = intersect(wt_c_of_max(l),wt_c_tuning_curve_of(l,:),'stable');
idxwt(l,1) = ia;
pref_ang_of(l,1) = hd_tuning_axis(idxwt(l));
end

figure()%%% WILDTYPE COMPRESSION SCATTER/HISTOPLOT
%Create x data histogram on top
g(1,1)=gramm('x',rad2deg(pref_ang_sq));
g(1,1).set_layout_options('Position',[0 0.8 0.8 0.2],... 
    'legend',false,... 
    'margin_height',[0.02 0.05],... 
    'margin_width',[0.1 0.02],...
    'redraw',false); 
g(1,1).set_names('x','','y','Neurons');
g(1,1).stat_bin('geom','line','fill','all','nbins',8); %histogram
g(1,1).axe_property('XTickLabel',''); 
%Create a scatter plot
g(2,1)=gramm('x',rad2deg(pref_ang_sq),'y',rad2deg(pref_ang_of),'color',wt_c_mean_fr_diff);
g(2,1).set_names('x','Preferred Angle in Compression (Degrees)','y','Preferred Angle in Baseline','color','');
g(2,1).geom_point(); 
g(2,1).set_point_options('base_size',10);
g(2,1).set_continuous_color('colormap','parula','Clim',[-8 8])
g(2,1).set_layout_options('Position',[0 0 0.8 0.8],...
    'legend',false,... 
    'margin_height',[0.1 0.02],...
    'margin_width',[0.1 0.02],...
    'redraw',false);
%Create y data histogram on the right
g(3,1)=gramm('x',rad2deg(pref_ang_of));
g(3,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
    'legend',false,...
    'margin_height',[0.1 0.02],...
    'margin_width',[0.02 0.05],...
    'redraw',false);
g(3,1).set_names('x','','y','Neurons');
g(3,1).stat_bin('geom','line','fill','all','nbins',8); %histogram
g(3,1).coord_flip();
g(3,1).axe_property('XTickLabel','');
g.set_color_options('map',[color.blue1])
g.draw();
colormap(brewermap([],'*RdYlBu'));

for j = 1:height(wt_e_of)
    wt_e_tuning_curve_of(j,:) = wt_e_of.Raw_Curve{j,1};
    wt_e_tuning_curve_sq(j,:) = wt_e_sq.Raw_Curve{j,1};
end
hd_tuning_axis = linspace(0,2*pi,60);

wt_e_sq_max = (max(wt_e_tuning_curve_sq,[],2));
wt_e_of_max = (max(wt_e_tuning_curve_of,[],2));

for l = 1:size(wt_e_sq_max,1)
[h,hg,ia] = intersect(wt_e_sq_max(l),wt_e_tuning_curve_sq(l,:),'stable');
idxwt(l,1) = ia;
pref_ange_sq(l,1) = hd_tuning_axis(idxwt(l));
[h,hg,ia] = intersect(wt_e_of_max(l),wt_e_tuning_curve_of(l,:),'stable');
idxwt(l,1) = ia;
pref_ange_of(l,1) = hd_tuning_axis(idxwt(l));
end
figure() %%% WILDTYPE EXPANSION SCATTER/HISTOPLOT
%Create x data histogram on top
g(1,1)=gramm('x',rad2deg(pref_ange_sq));
g(1,1).set_layout_options('Position',[0 0.8 0.8 0.2],... 
    'legend',false,... 
    'margin_height',[0.02 0.05],... 
    'margin_width',[0.1 0.02],...
    'redraw',false); 
g(1,1).set_names('x','','y','Neurons');
g(1,1).stat_bin('geom','line','fill','all','nbins',8); %histogram
g(1,1).axe_property('XTickLabel',''); 
%Create a scatter plot
g(2,1)=gramm('x',rad2deg(pref_ange_sq),'y',rad2deg(pref_ange_of),'color',wt_e_mean_fr_diff);
g(2,1).set_names('x','Preferred Angle in Expansion (Degrees)','y','Preferred Angle in Baseline','color','');
g(2,1).geom_point(); 
g(2,1).set_point_options('base_size',10);
g(2,1).axe_property('xlim',[0 360],'ylim',[0 360]);
g(2,1).set_continuous_color('colormap','parula','Clim',[-8 8])
g(2,1).set_layout_options('Position',[0 0 0.8 0.8],...
    'legend',false,... 
    'margin_height',[0.1 0.02],...
    'margin_width',[0.1 0.02],...
    'redraw',false);
%Create y data histogram on the right
g(3,1)=gramm('x',rad2deg(pref_ange_of));
g(3,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
    'legend',false,...
    'margin_height',[0.1 0.02],...
    'margin_width',[0.02 0.05],...
    'redraw',false);
g(3,1).set_names('x','','y','Neurons');
g(3,1).stat_bin('geom','line','fill','all','nbins',8); %histogram
g(3,1).coord_flip();
g(3,1).axe_property('XTickLabel','');
g.set_color_options('map',[color.orange1]);
g.draw();
colormap(brewermap([],'*RdYlBu'));