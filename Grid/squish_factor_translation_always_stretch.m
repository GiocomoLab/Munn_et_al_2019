function [rho_max,mean_rho, yShift_cm, xShift_cm, lambda] = squish_factor_translation_always_stretch(file1, file2, file3, file4,boxSize)
% Caitlin Mallory
% 1/12/18

% This script is used to compute the correlation between the grid map in
% the compressed environment and the grid map in the original, square
% environment. To determine how much "squishing" or grid deformation occured,
% the "squished" map is stretched in the direction of the compression in
% increasing steps between the width of the squished box, and the width of the
% original box. For each successive stretch of the squish map, the script computes the
% correlation between the resized squish map and the overlapping portion of
% the original map. In addition, for each resized squish map, the correlation is computed
% when shifting the X and Y position of the stretched-squish map in
% increments (the magnitude of shift allowed can be set by the user below).
% This is to allow for translation of the grid pattern.

% file1 = Spike file, open field
% file2 = Position file, open field
% file3 = Spike file, squish
% file4 = Position file, squish

%NOTES ON INTERPRETATION:
%if y is negative, the squish map has shifted "down" relative to the
%original.
%if y is positive, the squish map has shifted "up" relative to the
%original.
% if x is negative, the squish map has shifted "left" relative to the
% original
%if x is positive, the squish map has shifted "right" relative to the
%original

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters to set:
x_maxShift = 13; %bins (*2 to get cm)
y_maxShift = 13; %bins
minAreaForComparing = 0.25; % This is the minimum proportion of the original
%box size that must overlap with the stretched, translated squish map. Usually 0.35
%The rho values for translations that result in overlap lower than this
%value will be set to NaN and not considered.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rotate = 0; %start out assuming that no rotation is needed
plotfig = 0;
[map_rectangle, rotate] = map_maker_squish(file3,file4,1,rotate);
%run map maker on the squish session first to determine if a rotation
%occured. if so, will need to also rotate the baseline map!
[map_square, rotate] = map_maker_squish(file1,file2,0,rotate);


% find maps of squish and non-squish sessions
%map_no_squish = map_maker(file1,file2,0);
%map_squish = map_maker(file3,file4,1);

% find dimensions of maps
[x1,y1] = size(map_rectangle);
[x2,y2] = size(map_square);

% check to make sure number of rows is the same (sanity check)
if x1 ~= x2
    map_rectangle = map_rectangle(1:x2,1:y1);
end

rho = zeros(2*x_maxShift+1,2*y_maxShift+1,(y2-y1+1));
for x = -x_maxShift:1:x_maxShift
    %start by shifting the squish map left or right relative to the
    %original map.
    
    for y = -y_maxShift:1:y_maxShift
        % also shift the map in the up or down direction relative to the
        % original map
        
        for s = y1:y2
            
            map_resize = imresize(map_rectangle, [size(map_rectangle,1), s]); %resize squished map to be same number of rows as original, and new number of columns
            
            if x >= 0 && y >= 0
                map_rectangle_to_compare = map_resize(1:size(map_resize,1)-y, 1:min(size(map_square,2)-x,size(map_resize,2)));
                baseline_comparison = map_square(y+1:end,x+1:min(x + size(map_resize,2),size(map_square,2)));
                
            elseif x <= 0 && y >= 0
                map_rectangle_to_compare = map_resize(1:size(map_resize,1)-y,abs(x)+1:size(map_resize,2));
                baseline_comparison = map_square(y+1:end,1:min(size(map_resize,2)- abs(x),size(map_square,2)));
                
            elseif x >= 0 && y <= 0
                
                map_rectangle_to_compare = map_resize(abs(y)+1:end,1:min(size(map_square,2) - x,size(map_resize,2)));
                baseline_comparison = map_square(1:end-abs(y),x+1:min(size(map_square,2),x+size(map_resize,2)));
                
            elseif x <= 0 && y <= 0
                
                map_rectangle_to_compare = map_resize(abs(y)+1:end,abs(x)+1:end);
                baseline_comparison = map_square(1:size(map_square,1)-abs(y),1:size(map_resize,2)-abs(x));
                
            end
            
            %calculate correlation coefficient
            X = corr2(map_rectangle_to_compare,baseline_comparison);
            index1 = y + y_maxShift + 1;
            index2 = x + x_maxShift + 1;
            index3 = s-y1+1;
            
            % restrict the correlations so that the area being compared is
            % at least 1/3 the area of the original box.
            if size(baseline_comparison,1)*size(baseline_comparison,2) >= minAreaForComparing*size(map_square,1)*size(map_square,2)
                rho(index1,index2,index3) = X;
            else
                rho(index1,index2,index3) = nan;
            end
        end
    end
end

xBinWidth = 2;
yBinWidth  = 2;
[rho_max, linearIndexesOfMaxes] = max(rho(:));
[yShiftIndex,xShiftIndex,stretchIndex] = ind2sub(size(rho),linearIndexesOfMaxes);
lambda = (stretchIndex-1)/(y2-y1); %normalized squish factor
stretch_units = stretchIndex - 1;
mean_rho = nanmean(rho(:));
yShift = -y_maxShift + yShiftIndex - 1; %in bins
xShift = -x_maxShift + xShiftIndex - 1; %in bins
xShift_cm = xBinWidth*xShift;
yShift_cm = yBinWidth*yShift; %in cm
displambda = sprintf('%s%s','Stretch Factor (Lambda) = ',num2str(lambda));
dispyShift_cm = sprintf('%s%s','Y Shift (cm) = ',num2str(yShift_cm));
dispxShift_cm = sprintf('%s%s','X Shift (cm) = ',num2str(xShift_cm));
disprho_max = sprintf('%s%s','Correlation (Rho) = ',num2str(rho_max));
dispmean_rho = sprintf('%s%s','Correlation (Rho) = ',num2str(mean_rho));

disp(displambda);
disp(dispyShift_cm);
disp(dispxShift_cm);
disp(disprho_max);
end