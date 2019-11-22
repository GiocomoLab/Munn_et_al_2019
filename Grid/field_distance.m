function [dist] = field_distance(fieldPos)

N = size(fieldPos,1);
x = fieldPos(:,1);
y = fieldPos(:,2);

% Find the field closest to the centre
dist = sqrt(x.^2 + y.^2);
[minDist, minInd] = min(dist);
cField = [x(minInd), y(minInd)];

% Find distance between centre field and the rest of the fields
dist = zeros(N,1);
for ii = 1:N
   if ii ~= minInd
       dist(ii) = sqrt((x(ii)-cField(1))^2 + (y(ii)-cField(2))^2);
   else
       % Set the distance to the centre field itself to infinity
       dist(ii) = inf;
   end
end