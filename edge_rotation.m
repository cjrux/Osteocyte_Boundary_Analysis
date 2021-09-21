function [border,major,minor,c,orientation,ang,Xc,Yc]= edge_rotation(border,s,I,mod1)

%hold on
major=s(I).MajorAxisLength; %determining major axis length
minor=s(I).MinorAxisLength; %determining minor axis length
c=s(I).Centroid; %determining centroid position
orientation=s(I).Orientation; %determining angle of "ellipse" from x-axis
%scatter(c(1,1),c(1,2),'g','*','LineWidth', 2); %plotting centroid
%hold on
% Create rotated outlines of lacuna boundary to prep and perform dilation
%hold on
%scatter(c(1,1),c(1,2),'g','*','LineWidth', 2);
%hold on
orientation_slope=sind(orientation)/cosd(orientation); %converting orientation into slope
b=c(1,2)-orientation_slope*c(1,1);%finding y intercept of orientation line
norm_slope=-inv(orientation_slope);%finding slope normal to orientation line

X = border(:,1);
Y = border(:,2);
ang = 270-orientation;
% Specify the coordinates of the center of rotation
Xc = c(1,1) ;  % Rotate about the 1/4 chord point
Yc = c(1,2) ;
% The data is roated in a three-step process
% Step 1) Shift the data to the rotation center
Xs = X - Xc;  % shifted data
Ys = Y - Yc;
% Step 2) Rotate the data
Xsr =  Xs*cosd(ang) + Ys*sind(ang);    % shifted and rotated data
Ysr = -Xs*sind(ang) + Ys*cosd(ang);    %
% Step 3) Un-shift the data (back to the original coordinate system)
Xr = Xsr + Xc;  % Rotated data
Yr = Ysr + Yc;

%For verification of boundary, can plot
%plot(Xr,Yr,'m','Linewidth',2)
% xlim([0 length(mod1)])
% ylim([0 length(mod1)])
%newborder(:,1) =  border(:,1)*cosd(orientation) + border(:,2)*sind(orientation);    % shifted and rotated data
%newborder(:,2) = -border(:,1)*sind(orientation) + border(:,2)*cosd(orientation);
border(:,1)=Xr;
border(:,2)=Yr;
hold on