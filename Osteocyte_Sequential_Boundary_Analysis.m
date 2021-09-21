clear all; close all; clc;
filename='test_file.csv';%input filename for .csv file
mod = csvread(filename,2,0);%reading correct columns and rows for data 2nd row, and all columns -- may need to change
%% Masking Stepm -- setting all points within lacunae equal to zero
element_size=5;%erosion element size
lower_bound=0;%data lower bound
upper_bound=30;%data upperbound
[modx,inner]=mask(mod,element_size,lower_bound,upper_bound);
area_points=length(modx)^2-sum(inner,'all');
%% Actual Lacunar Border -- points used for actual appropriately eroded lacuna edge
element_size=7; %best to capture whole osteocyte, typically 4 -- determines intensity of erosion
[border2,mod2,s2,I2]=edgefinder(modx,element_size,0,30); %image processing function
hold on
scatter(border2(:,1),border2(:,2),'c.') %can check border overlayed on mod map here
%% Over-Eroded Border -- points used for dilation and convolution matrix
element_size=10; %10 good value at 0.5 um dilations -- determines intensity of erosion
[border,mod1,s,I]=edgefinder(modx,element_size,0,30); %image processing function
hold on
%plot(border(:,1),border(:,2),'r','Linewidth',2)
scatter(border(:,1),border(:,2),'c.','Linewidth',2) %can check border overlayed on mod map here
%% Rotating lacunar edge and map for later use to create non-distorted dilations
[border,major,minor,c,orientation,ang,Xc,Yc]= edge_rotation(border,s,I,mod1);
[border2,major2,minor2,c2,orientation2,ang2,Xc2,Yc2]= edge_rotation(border2,s2,I2,mod1);
[x_rotated,y_rotated,modt]=rotate(mod1,ang2,Xc2,Yc2); %function for rotating modulus matrix
%% --Dilation and Creation of Sequential Boundaries--
%using dilation function to create dilated boundaries for 'binning' data
num_dilations=10; %must choose even number...could change this if wanted
map_dim=12; %um
step_size=0.2; %um
boundmat = dilation(border, major, minor,c,num_dilations,mod1,map_dim,step_size); %dilation founding
new_bound=boundmat{1}; %assigning first boundary to matrix for later use
%scatter(new_bound(:,1),new_bound(:,2),'.','k','Linewidth',2)
%% --Visualize Map Rotation and Unsmoothed Sequential Boundaries--
%Can show figure if interested in rotation and dilation fit
figure(2)
scatter(x_rotated, y_rotated, 5, modt);
colormap(hot)
colorbar;
    a=30; %colorbar upper bound
    b=0; %colorbar lower bound
    caxis([b a])%these values can be changed depending on modulus range
xlim([0 length(mod1)])
ylim([0 length(mod1)])
axis equal
hold on
scatter(border(:,1),border(:,2),'.','b','Linewidth',8)
hold on
T=cell2mat(boundmat)
X_border = T(:, 1:2:end);
Y_border = T(:, 2:2:end);
scatter(X_border(:),Y_border(:),'.','b','Linewidth',10)
%plot(X_border(:),Y_border(:),'.','b','Linewidth',10)
%% --Smoothing Boundary Lines--
start_size=11; %original is 21 -- For large initial smooth adjust this value -- this is kernel size
percent_increase=1.05; %original is 1.2 To speed up smoothing from one boundary to another increase this value
finalmat = smooth(boundmat,mod1,start_size,percent_increase); %smoothing function
%% --Putting Points in Counterclockwise Order With Cell Array Function--
[x2,y2]= counterclockwise_smooth(Xc, Yc, border2,finalmat); %putting adjacent boundaries in CCW and CW orders for inpolygon use
% Visualize New Boundary -- Not necessary
scatter(x_rotated, y_rotated, 5, modt);
colormap(hot)
colorbar;
    a=30; %colorbar upper bound
    b=0; %colorbar lower bound
    caxis([b a])%these values can be changed depending on modulus range
xlim([0 length(mod1)])
ylim([0 length(mod1)])
hold on
%graphing smoothed boundaries on modulus map
for i=1:length(x2)
    pointx=x2{i}
    pointy=y2{i}
    scatter(pointx,pointy,'b.')
    hold on
end
hold on
scatter(border2(:,1),border2(:,2),'b.')
%% -- Prepping for inpolygon and cataloging points within given boundaries
%creating matrix columns and cataloging points within two concentric boundaries 
 for j=1:length(x2)-1
    x=x2{1,j};
    y=y2{1,j};
    x1=x2{1,j+1};
    y1=y2{1,j+1};
 xv = [rot90(x,3) NaN rot90(x1)];
 yv = [rot90(y,3) NaN rot90(y1)];
in = inpolygon(x_rotated,y_rotated,xv,yv);
%hold on
%plot(x_rotated(in),y_rotated(in),'g.') % points inside
%x_in{1,j}=x_rotated(in);
%y_in{1,j}=y_rotated(in);
points_in{j}=in;
 end
%% --Concatenate points within polygon into array
for i=1:length(points_in)
hold=points_in{i};
x_holder=x_rotated(hold);
y_holder=y_rotated(hold);
x_points{i}=x_holder;
y_points{i}=y_holder;

modmapx{i}=x_rotated;
modmapy{i}=y_rotated;
end
%% -- Matching points with given modulus at point and prepping for histograms--
% Determining which points are of interest
for i=1:length(x_points)
Liax = ismember(modmapx{i},x_points{i});
Liay = ismember(modmapy{i},y_points{i});
Loc=Liax.*Liay;
E=find(Loc);
Hmod{i}=modt(E);
%could add step to filter data under certain threshold i.e. 5 GPa
% indices = find(Hmod>0);
% Hmod_new=Hmod(indices);
vector = Hmod{i};
nvector = vector(vector>5 & vector<90);%can change threshold
r(i)=range(Hmod{i});%determining range
pd = fitdist(nvector,'Normal')%fitting to normal distribution
maximums(i)=ceil(max(Hmod{i}));%determining max value and rounding up to next whole integer value
%creating table of relevant statistical properties
stat_values(i,1)=range(nvector);
stat_values(i,2)=median(nvector);
stat_values(i,3)=mean(nvector);
stat_values(i,4)=std(nvector);
stat_values(i,5)=2*sqrt(2*pd.sigma^2*log(2));
num_points(i,1)=length(nvector);
Hmod{i}=nvector;
% stat_values(i,1)=range(Hmod{i});
% stat_values(i,2)=median(Hmod{i});
% stat_values(i,3)=mean(Hmod{i});
% stat_values(i,4)=std(Hmod{i});
% stat_values(i,5)=2*sqrt(2*pd.sigma^2*log(2))
clear Loc
clear E
clear Liax
clear Liay
end
colNames = {'range','median','mean','standard deviation','FWHM'};
sTable = array2table(stat_values,'VariableNames',colNames)

%Can check binning is correct with scatter
%hold on
%scatter(x_points{2},y_points{2},'.c')

%% Creating Histograms
figure (4)
for i=1:length(Hmod);
subplot(2,length(Hmod)/2,i)
bin_edges=(0:maximums(i));
%bin_edges=(1:maximums(i));
h=histogram(Hmod{i},bin_edges)
xlim([0, maximums(i)]);
%xlim([0, maximums(i)]);
[maxcount, whichbin] = max(h.Values);
title("Maximum occurs in bin " + whichbin + ...
      " between " + h.BinEdges(whichbin) + " and " + h.BinEdges(whichbin+1));
xlabel('Moduli')
%title('PM02 0.25 um dilation  # bins= discrete range')
end
sTable %printing table of statistics for each concentric region
