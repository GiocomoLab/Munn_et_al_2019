function [Results] = direction_code_v2(sesh,sesh_sq,pos,pos_sq,animal,plot_figure,type)
Results = table();
% plot_figure = 1; %% IF YOU WANT TO PLOT AND SAVE THE TUNING CURVES SET 1, OTHERWISE SET 0
for n = 1:length(sesh)
    
    spikefile_of = sesh{n};
    posfile_of = pos{n};
    posfile_sq = pos_sq{n};
    spikefile_sq = sesh_sq{n};
    
    
    load(posfile_sq)
    load(spikefile_sq)
    
    posy_sq = -posy;
    posy2_sq = -posy2; %the output from database maker is the mirror image of tint. This flips the output to be the same
    posx_sq = posx;
    posx2_sq = posx2;
    post_sq = post;
    cellTS_sq = cellTS;
    minX = nanmin(posx_sq); maxX = nanmax(posx_sq);
    minY = nanmin(posy_sq); maxY = nanmax(posy_sq);
    xLength = maxX - minX;
    yLength = maxY - minY;
    posx_sq = (posx_sq - min(posx_sq)); % Scale to boxsize (100cm)
    posy_sq = (posy_sq - min(posy_sq));
    posx2_sq = (posx2_sq - min(posx2_sq)); % Scale to boxsize (100cm)
    posy2_sq = (posy2_sq - min(posy2_sq));
    if xLength > yLength
        scalefac = 100/xLength;
        posx_sq = posx_sq*scalefac;
        posx2_sq = posx2_sq*scalefac;
        posy_sq = posy_sq*scalefac;
        posy2_sq = posy2_sq*scalefac;
    elseif yLength > xLength
        scalefac = 100/yLength;
        posx_sq = posx_sq*scalefac;
        posx2_sq = posx2_sq*scalefac;
        posy_sq = posy_sq*scalefac;
        posy2_sq = posy2_sq*scalefac;
    end
    minX = nanmin(posx_sq); maxX = nanmax(posx_sq);
    minY = nanmin(posy_sq); maxY = nanmax(posy_sq);
    xLength = maxX - minX;
    yLength = maxY - minY;
    if yLength > xLength+10 % if the stable wall (cue wall) is west % rotate 90 degrees ccw
        posx_new = -posy_sq; posx2_new = -posy2_sq; posy_sq = posx_sq; posy2_sq = posx2_sq;
        posx_sq = posx_new; posx2_sq = posx2_new;
        if min(posx_sq) < 0
            posx_sq = posx_sq + abs(min(posx_sq));
            posx2_sq = posx2_sq+ abs(min(posx2_sq));
        end
        rotation_sq = 1;
    else
        rotation_sq = 0;
    end
    clear posx posy post cellTS
    load(posfile_of)
    load(spikefile_of)
    posy_of = -posy;
    posy2_of = -posy2; %the output from database maker is the mirror image of tint. This flips the output to be the same
    posx_of = posx;
    posx2_of = posx2;
    post_of = post;
    cellTS_of = cellTS;
    minX = nanmin(posx_of); maxX = nanmax(posx_of);
    minY = nanmin(posy_of); maxY = nanmax(posy_of);
    xLength = maxX - minX;
    yLength = maxY - minY;
    posx_of = (posx_of - min(posx_of)); % Scale to boxsize (100cm)
    posy_of = (posy_of - min(posy_of));
    posx2_of = (posx2_of - min(posx2_of)); % Scale to boxsize (100cm)
    posy2_of = (posy2_of - min(posy2_of));
    if xLength > yLength
        scalefac = 100/xLength;
        posx_of = posx_of*scalefac;
        posx2_of = posx2_of*scalefac;
        posy_of = posy_of*scalefac;
        posy2_of = posy2_of*scalefac;
    elseif yLength > xLength
        scalefac = 100/yLength;
        posx_of = posx_of*scalefac;
        posx2_of = posx2_of*scalefac;
        posy_of = posy_of*scalefac;
        posy2_of = posy2_of*scalefac;
    end
    minX = nanmin(posx_of); maxX = nanmax(posx_of);
    minY = nanmin(posy_of); maxY = nanmax(posy_of);
    xLength = maxX - minX;
    yLength = maxY - minY;
    
    if rotation_sq == 1 % if the stable wall (cue wall) is west % rotate 90 degrees ccw
        posx_new = -posy_of; posx2_new = -posy2_of; posy_of = posx_of; posy2_of = posx2_of;
        posx_of = posx_new; posx2_of = posx2_new;
        if min(posx_of) < 0
            posx_of = posx_of + abs(min(posx_of));
            posx2_of = posx2_of + abs(min(posx2_of));
        end
        rotation_of = 1;
    else
        rotation_of = 0;
    end
    
    if rotation_sq == 1 && rotation_of == 1
        disp('Rotated both the maps')
    elseif rotation_sq == 0
        disp('No Squish Rotation')
    elseif rotation_sq == 1 && rotation_of == 0
        disp('Something went wrong with the squish rotations')
    end
    if type == 2
        if max(posy_of) > max(posx_of)
            keyboard
        end
    end
    lowSpeedThreshold = 2;
    highSpeedThreshold = 100;
    
    speed_of = speed2D(posx_of,posy_of,post_of);
    speed_sq = speed2D(posx_sq,posy_sq,post_sq);
    bad_ind_of = find(speed_of < lowSpeedThreshold | speed_of > highSpeedThreshold);
    bad_ind_sq = find(speed_sq < lowSpeedThreshold | speed_sq > highSpeedThreshold);
    
    posx_of(bad_ind_of) = NaN;
    posy_of(bad_ind_of) = NaN;
    posx_sq(bad_ind_sq) = NaN;
    posy_sq(bad_ind_sq) = NaN;
    speed_of(bad_ind_of) = NaN;
    speed_sq(bad_ind_sq) = NaN;
    
     if length(posy2_sq) > length(posy_sq)
                posy2_sq = posy2_sq(1:length(posy_sq));
     end
      if length(posy2) > length(posy)
                posy2 = posy2(1:length(posy));
      end
        if length(posx2_sq) > length(posx_sq)
                posx2_sq = posx2_sq(1:length(posx_sq));
     end
      if length(posx2) > length(posx)
                posx2 = posx2(1:length(posx));
     end
     
    roundCellTS = round(cellTS*100)/100;
    roundCellTS = round(cellTS_sq*100)/100;
    ind_post = round(roundCellTS*50)+1; %assumes 20 ms sampling frequency
    ind_post_sq = round(roundCellTS*50)+1;
    ind_post(ind_post>numel(post)) = [];
    ind_post_sq(ind_post_sq>numel(post_sq)) = [];
    
    % take out spikes that occur during stationary points
    ind_post_of = ind_post(~ismember(ind_post,bad_ind_of));
    ind_post_sq = ind_post_sq(~ismember(ind_post_sq,bad_ind_sq));
    
    % find fr as a function of heading direction
    try
    [hdAxis , hd_fr] = computeMVL(posx_of,posx2_of,posy_of,posy2_of,ind_post_of);
    
    [hdAxis_sq , hd_fr_sq] = computeMVL(posx_sq,posx2_sq,posy_sq,posy2_sq,ind_post_sq);
    catch
        keyboard
     end
    
        
    % compute mean vector length and argument
    exp_dir = exp(-1i*hdAxis);
    rayleigh_vector = pi/(numel(hdAxis)*sin(pi/numel(hdAxis)))*sum(hd_fr.*exp_dir)/sum(hd_fr);
    mvl = sqrt(real(rayleigh_vector)^2+imag(rayleigh_vector)^2);
    mv_arg = -atan2(imag(rayleigh_vector),real(rayleigh_vector))*180/pi;
    
    if mv_arg < 0
        mv_arg = mv_arg+360; %make mv_arg give proper degree output
    else
    end
    
    exp_dir_sq = exp(-1i*hdAxis_sq);
    rayleigh_vector_sq = pi/(numel(hdAxis_sq)*sin(pi/numel(hdAxis_sq)))*sum(hd_fr_sq.*exp_dir_sq)/sum(hd_fr_sq);
    mvl_sq = sqrt(real(rayleigh_vector_sq)^2+imag(rayleigh_vector_sq)^2);
    mv_arg_sq = -atan2(imag(rayleigh_vector_sq),real(rayleigh_vector_sq))*180/pi;
    if mv_arg_sq < 0
        mv_arg_sq = mv_arg_sq+360; %make mv_arg give proper degree output
    else
    end
    % compute peak firing rate and preferred angle
    [peak_rate,pref_angle] = nanmax(hd_fr);
    [peak_rate_sq,pref_angle_sq] = nanmax(hd_fr_sq);
    
    % compute the half-width, half-mean
    half_mean = peak_rate/2;
    [~,half_inds_left] = min(abs(half_mean - hd_fr(1:pref_angle)));
    [~,half_inds_right] = min(abs(half_mean - hd_fr(pref_angle:end)));
    hw_hm = (half_inds_right +pref_angle - 1 - half_inds_left)*(360/numel(hd_fr)); % 6 degrees per bin
    pref_angle_radians = (hdAxis(pref_angle));
    
    half_mean_sq = peak_rate_sq/2;
    [~,half_inds_left_sq] = min(abs(half_mean_sq - hd_fr_sq(1:pref_angle_sq)));
    [~,half_inds_right_sq] = min(abs(half_mean_sq - hd_fr_sq(pref_angle_sq:end)));
    hw_hm_sq = (half_inds_right_sq +pref_angle_sq - 1 - half_inds_left_sq)*(360/numel(hd_fr_sq)); % 6 degrees per bin
    pref_angle_radians_sq = (hdAxis_sq(pref_angle_sq));
    
    % compute the stability
    % divide the session into quarters
    % numpoints = numel(posx_of); div_axis = round(linspace(1,numpoints,5));
    % hd_fr_mat = nan(4,numel(hdAxis));
    % for k = 1:4
    %     spk_ind = ind_post(ind_post >= div_axis(k) & ind_post <= div_axis(k+1))  - div_axis(k) + 1 ;
    %     [~,hd_fr_mat(k,:)] = computeMVL(posx_of(div_axis(k):div_axis(k+1)), ...
    %     posx2_of(div_axis(k):div_axis(k+1)),posy_of(div_axis(k):div_axis(k+1)),...
    %     posy2_of(div_axis(k):div_axis(k+1)),spk_ind);
    % end
    % keyboard
    % numpoints_sq = numel(posx_sq); div_axis_sq = round(linspace(1,numpoints_sq,5));
    % hd_fr_mat_sq = nan(4,numel(hdAxis_sq));
    % for k = 1:4
    %     spk_ind_sq = ind_post_sq(ind_post_sq >= div_axis_sq(k) & ind_post_sq <= div_axis_sq(k+1))  - div_axis(k) + 1 ;
    %     [~,hd_fr_mat_sq(k,:)] = computeMVL(posx_sq(div_axis_sq(k):div_axis_sq(k+1)), ...
    %     posx2_sq(div_axis_sq(k):div_axis_sq(k+1)),posy(div_axis_sq(k):div_axis_sq(k+1)),...
    %     posy2_sq(div_axis_sq(k):div_axis_sq(k+1)),spk_ind_sq);
    % end
    %
    % % compare all
    % combinations = [ 1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    % corr_all = nan(length(combinations),1);
    % for k = 1:length(combinations)
    %     compare = combinations(k,:);
    %     first = (hd_fr_mat(compare(1),:));
    %     second = (hd_fr_mat(compare(2),:));
    %     idx = (~isnan(second) & ~isnan(first));
    %     second = second(idx);
    %     first = first(idx);
    %     corr_all(k) = corr(first',second');
    % end
    % Stability = mean(corr_all);
    %
    % if isempty(Stability)
    %     keyboard
    %     Stability = NaN;
    % end
    %
    % combinations = [ 1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    % corr_all_sq = nan(length(combinations),1);
    % for k = 1:length(combinations)
    %     compare = combinations(k,:);
    %     first_sq = (hd_fr_mat_sq(compare(1),:));
    %     second_sq = (hd_fr_mat_sq(compare(2),:));
    %     idx_sq = (~isnan(second_sq) & ~isnan(first_sq));
    %     second_sq = second_sq(idx_sq);
    %     first_sq = first_sq(idx_sq);
    %     corr_all_sq(k) = corr(first_sq',second_sq');
    % end
    % Stability_sq = mean(corr_all_sq);
    %
    % if isempty(Stability_sq)
    %     keyboard
    %     Stability_sq = NaN;
    % end

    if iscell(animal)
    Results.animal{n,1} = animal{n};
    else
        Results.animal{n,1} = sprintf('%s',animal);
    end
    Results.session{n,1} = sesh{n};
    Results.mvl(n,1) = mvl;
    Results.mv_arg(n,1) = mv_arg;
    Results.pref_angle(n,1) = pref_angle_radians;
    Results.peakFR(n,1) = peak_rate;
    %  Results.stability(n) = Stability;
    Results.Raw_Curve{n,1} = hd_fr;
    Results.half_width(n,1) = hw_hm;
    if iscell(animal)
    Results.animal_sq{n,1} = animal{n};
    else
        Results.animal_sq{n,1} = sprintf('%s',animal);
    end
    Results.session_sq{n,1} = sesh_sq{n};
    Results.mvl_sq(n,1) = mvl_sq;
    Results.mv_arg_sq(n,1) = mv_arg_sq;
    Results.pref_angle_sq(n,1) = pref_angle_radians_sq;
    Results.peakFR_sq(n,1) = peak_rate_sq;
    %  Results.stability_sq(n) = Stability_sq;
    Results.Raw_Curve_sq{n,1} = hd_fr_sq;
    Results.half_width_sq(n,1) = hw_hm_sq;
    
    hdDir_of = calcHeadDirection(posx_of,posy_of,posx2_of,posy2_of);
    hdDir_sq = calcHeadDirection(posx_sq,posy_sq,posx2_sq,posy2_sq);
    % Find the direction of the rat at the spike times
    spkDir_of = hdDir_of(ind_post_of);
    spkDir_of(isnan(spkDir_of)) = [];
    spkDir_sq = hdDir_sq(ind_post_sq);
    spkDir_sq(isnan(spkDir_sq)) = [];
    p.hdBinWidth = 6; % [degrees]
    p.hdTimeBinWidth = 6; % [degrees]
    p.hdSmoothingMode = 1;
    p.hdSmoothingWindowSize = 14.5; % [degrees]
    p.percentile = 50; % [%]
    p.hdAlphaValue = 10000;
    % Do the head direction analysis
    hd_of = hdstat(spkDir_of*2*pi/360,hdDir_of*2*pi/360,0.02, 1);
    [hdMapA_of, dirPDF_of] = hdMapAdaptiveSmoothing(spkDir_of, hdDir_of, 0.02, p);
    hd_sq = hdstat(spkDir_sq*2*pi/360,hdDir_sq*2*pi/360,0.02, 1);
    [hdMapA_sq, dirPDF_sq] = hdMapAdaptiveSmoothing(spkDir_sq, hdDir_sq, 0.02, p);
    % plot head direction map
    l = length(hdMapA_of);
    angles = 0:2*pi/l:2*pi;
    angles = angles';
    map = hdMapA_of;
    map(end+1) = hdMapA_of(1);
    % X and Y coordinates of the head direction rate map
    plotRX = cos(angles) .* map;
    plotRY = sin(angles) .* map;
    minX = min(plotRX);
    maxX = max(plotRX);
    minY = min(plotRY);
    maxY = max(plotRY);
    minValue = min([minX, minY]);
    maxValue = max([maxX, maxY]);
    maxValue = max([abs(minValue), abs(maxValue)]);
    % Add the axis
    aLength = 2* maxValue;
    addL = 0.05 * aLength;
    if plot_figure
        figure(11)
        cla
        hold off
        line([-maxValue - addL, maxValue + addL], [0,0], 'color', [0.5, 0.5, 0.5]);
        line([0,0] ,[-maxValue - addL, maxValue + addL], 'color', [0.5, 0.5, 0.5]);
        hold on;
        % Plot the rate map
        plot(plotRX, plotRY, 'k','lineWidth',1);
        hold off
        axis off
        axis image
        recording = strsplit(spikefile_of,'\');   
        recording  = recording{end};
        fName = sprintf('%s%s%s',animal{n},'_',recording,'_HeadDirectionMapAdaptiveSmoothed','.jpg');
        saveas(gcf,fName,'jpg');
    end
    clear map plotRX plotRY
    %%
    l = length(hdMapA_sq);
    angles = 0:2*pi/l:2*pi;
    angles = angles';
    map = hdMapA_sq;
    map(end+1) = hdMapA_sq(1);
    % X and Y coordinates of the head direction rate map
    plotRX = cos(angles) .* map;
    plotRY = sin(angles) .* map;
    minX = min(plotRX);
    maxX = max(plotRX);
    minY = min(plotRY);
    maxY = max(plotRY);
    minValue = min([minX, minY]);
    maxValue = max([maxX, maxY]);
    maxValue = max([abs(minValue), abs(maxValue)]);
    % Add the axis
    aLength = 2* maxValue;
    addL = 0.05 * aLength;
    if plot_figure
        figure(11)
        cla
        hold off
        line([-maxValue - addL, maxValue + addL], [0,0], 'color', [0.5, 0.5, 0.5]);
        line([0,0] ,[-maxValue - addL, maxValue + addL], 'color', [0.5, 0.5, 0.5]);
        hold on;
        % Plot the rate map
        plot(plotRX, plotRY, 'k','lineWidth',1);
        hold off
        axis off
        axis image
        recording = strsplit(spikefile_sq,'\');
        recording  = recording{end};
        fName = sprintf('%s%s%s',animal{n},'_',recording,'_HeadDirectionMapAdaptiveSmoothed','.jpg');
        saveas(gcf,fName,'jpg');
    end
    %%
    if plot_figure
        figure();
        recording = strsplit(spikefile_of,'\');
        recording  = recording{end};
        plotname = sprintf('%s%s',animal{n},'_',recording,'.jpg');
        figname = sprintf('%s%s',animal{n},'_',recording,'.fig');
        s = polarplot(hdAxis,hd_fr);
        saveas(gcf,plotname,'jpg');
        savefig(gcf,figname);
        close all
        figure();
        recording_sq = strsplit(spikefile_sq,'\');
        recording_sq = recording_sq{end};
        plotname_sq = sprintf('%s%s',animal{n},'_',recording_sq,'.jpg');
        figname_sq = sprintf('%s%s',animal{n},'_',recording_sq,'.fig');
        sq = polarplot(hdAxis_sq,hd_fr_sq);
        saveas(gcf,plotname_sq,'jpg');
        savefig(gcf,figname_sq);
        close all
    end
end
end
