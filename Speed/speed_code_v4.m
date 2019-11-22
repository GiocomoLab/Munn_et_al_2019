function [Results] = speed_code_v5(spikefile,posfile,animal,plot_figure,type)
Results = table();
% plot_figure = 1; %% IF YOU WANT TO PLOT AND SAVE THE TUNING CURVES SET 1, OTHERWISE SET 0
for n = 1:length(sesh)
    
    spikefile_of = sesh{n};
    posfile_of = pos{n};
    posfile_sq = pos_sq{n};
    spikefile_sq = sesh_sq{n};
 
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
    
    if yLength > xLength+10 % if the stable wall (cue wall) is west % rotate 90 degrees ccw
        posx_new = -posy_of; posx2_new = -posy2_of; posy_of = posx_of; posy2_of = posx2_sq;
        posx_sq = posx_new; posx2_sq = posx2_new;
        if min(posx_sq) < 0
            posx_sq = posx_sq + abs(min(posx_sq));
            posx2_sq = posx2_sq+ abs(min(posx2_sq));
        end
        rotation_sq = 1;
    else
        rotation_sq = 0;
    end
    
   
    lowSpeedThreshold = 2;
    highSpeedThreshold = 100;
    
    speed = speed2D(posx,posy,post);
  
    bad_ind = find(speed < lowSpeedThreshold | speed > highSpeedThreshold);
  
    
    posx(bad_ind) = NaN;
    posy(bad_ind) = NaN;
    speed(bad_ind) = NaN;

    
    
%compute speed at every time point
velx = diff([posx(1); posx]); vely = diff([posy(1); posy]); dt = 0.02;
speed = sqrt(velx.^2+vely.^2)/dt;
speedx = abs(velx/dt);
speedy = abs(vely/dt);


binWidth = 2;
speedVec = 5:binWidth:50;
sessionInd = round(linspace(1,numel(post),5));
sessionInd_sq = round(linspace(1,numel(post_sq),5));

[spiketrain] = computeFR(cellTS,post);
[spiketrain_sq] = computeFR(cellTS_sq,post_sq);
 %the tracking sampling rate is 50Hz and thus 20ms frames. The speed has the same                                                                            %sampling rate
fr = gauss_smoothing(spiketrain,20)*50;  % see attachment.
fr_sq = gauss_smoothing(spiketrain_sq,20)*50;
    
speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
speed = gauss_smoothing(speed,10); % gaussian kernel smoothing
speedx(isnan(speedx)) = interp1(find(~isnan(speedx)), speedx(~isnan(speedx)), find(isnan(speedx)), 'pchip'); % interpolate NaNs
speedx = gauss_smoothing(speedx,10); % gaussian kernel smoothing
speedy(isnan(speedy)) = interp1(find(~isnan(speedy)), speedy(~isnan(speedy)), find(isnan(speedy)), 'pchip'); % interpolate NaNs
speedy = gauss_smoothing(speedy,10); % gaussian kernel smoothing

meanFR = numel(cellTS)/max(post); % mean firing rate over whole session

% Bin speed for tuning curve
select = speed >= 2 & speed <= 100;
selectx = speedx >= 2 & speedx <= 100;
selecty = speedy >= 2 & speedy <= 100;
speedFilt = speed(select);
speedxFilt = speedx(selectx);
speedyFilt = speedy(selecty);
frFilt = fr(select);
frFiltx = fr(selectx);
frFilty = fr(selecty);

try
speedScore =corr(speedFilt,frFilt);
speedScore_x =corr(speedxFilt,frFiltx);
speedScore_y =corr(speedyFilt,frFilty);
catch
    keyboard
end

mean_speed = mean(speedFilt);
mean_speed_x = mean(speedxFilt);
mean_speed_y = mean(speedyFilt);

%% OF SLOPE AND INT
fone = polyfit(speedFilt,frFilt,1);
slope = fone(1);
intercept = fone(2);
fonex = polyfit(speedxFilt,frFiltx,1);
slopex = fonex(1);
interceptx = fonex(2);
foney = polyfit(speedyFilt,frFilty,1);
slopey = foney(1);
intercepty = foney(2);

%% OF STABILITY
fr_speed_maps = nan(4,numel(speedVec)-1);
for k = 1:4   
    speed_k = speed(sessionInd(k):sessionInd(k+1));
    fr_k = fr(sessionInd(k):sessionInd(k+1));
    
    for j = 1:numel(speedVec)-1
        start = speedVec(j); stop = speedVec(j+1);
        fr_speed_maps(k,j) = nanmean(fr_k(speed_k > start & speed_k < stop));
    end  
end
% take away any nan's
[~,badCol] = find(isnan(fr_speed_maps));
fr_speed_maps = fr_speed_maps(:,setdiff(1:numel(speedVec)-1,unique(badCol)));
% compute correlation for each comparison
session1 = [1 1 1 2 2 3]; session2 = [2 3 4 3 4 4];
correlations = nan(6,1);
for m = 1:6
    correlations(m) = corr(fr_speed_maps(session1(m),:)',fr_speed_maps(session2(m),:)');
end
stability  = nanmean(correlations);



Results.animal{n} = animal;
Results.session{n} = spikefile;
Results.score(n,1) = speedScore;
Results.scorex(n,1) = speedScore_x;
Results.scorey(n,1) = speedScore_y;
Results.slope(n,1) = slope;
Results.intercept(n,1) = intercept;
Results.slope_x(n,1) = slopex;
Results.intercept_x(n,1) = interceptx;
Results.slope_y(n,1) = slopey;
Results.intercept_y(n,1) = intercepty;
Results.mean_speed(n,1) = mean_speed;
Results.mean_x_speed(n,1) = mean_speed_x;
Results.mean_y_speed(n,1) = mean_speed_y;
Results.stability(n,1) = stability;

end
sesh
if plot_figure
    figure();
    recording = strsplit(spikefile,'\');
    recording  = recording{end};
    plotname = sprintf('%s%s',this_animal,'_',recording,'.fig');
    s = plot(speedAxis,speed_tuning);
    xlim([0 40]);
%     hold on
%     regline = refline(slope,intercept);
%     regline.Color = 'r';
%     regline.LineStyle = '--';
%     regline.LineWidth = 1.5;
    xlabel('Running Speed (cm/s)');
    ylabel('Firing Rate (Hz)');
    savefig(gcf,plotname);
    close all
end
end
