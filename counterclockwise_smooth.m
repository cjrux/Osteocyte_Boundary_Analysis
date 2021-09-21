function [x_cell,y_cell]= counterclockwise_smooth(Xc, Yc, border, finalmat)
% Putting points in CCW order using different method
xCenter = Xc;
yCenter = Yc;
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
%%
%----Outer boundary
for k=1:length(finalmat)
    mat=finalmat{1,k};
    X_border=mat(:,1);
    Y_border=mat(:,2);
for j=1:length(X_border)
x2=X_border(:,1);
y2=Y_border(:,1);
%then calculate the angles from the center
for i=1:length(x2)
angles2(i) = atan2d((y2(i)-yCenter),(x2(i)-xCenter));
end
%Then sort
[sortedAngles2, sortIndexes2] = sort(angles2);
x2_b = x2(sortIndexes2);  % Reorder x and y with the new sort order.
y2_b = y2(sortIndexes2);
%x2_b{1,i} = x2(sortIndexes2);  % Reorder x and y with the new sort order.
%y2_b{1,i} = y2(sortIndexes2);
clear angles2
end
x_cell{1}=x1;
y_cell{1}=y1;
x_cell{k+1}=x2_b;
y_cell{k+1}=y2_b;
end
end
%for i=1:length(x_cell)