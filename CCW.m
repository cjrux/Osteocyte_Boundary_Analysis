function [x1,y1]= CCW(centroid, border)
% Putting points in CCW order using different method
xCenter = centroid(:,1);
yCenter = centroid(:,2);
%---Inner boundary
x1=border(:,1);
y1=border(:,2);
%then calculate the angles from the center
for i=1:length(x1)
angles1(i) = atan2d((y1(i)-yCenter),(x1(i)-xCenter));
end
%Then sort
[sortedAngles1, sortIndexes1] = sort(angles1);
x1 = x1(sortIndexes1);  % Reorder x and y with the new sort order.
y1 = y1(sortIndexes1);
