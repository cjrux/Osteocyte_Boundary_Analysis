
function[x_rotated,y_rotated,modt]=rotate(mod1,ang,Xc,Yc)
n=length(mod1); %determining size of matrix 
xx=ones([1,n]);
xt=zeros;
for i=1:n
    xn=xx*i;
    xt=[xt, xn];
end
xt=rot90(xt,3);
xt=xt(2:end);
yt=zeros;
yy=1:n;
for i=1:n
    yt=[yt, yy];
    
end
yt=rot90(yt,3);
yt=yt(2:end);
% Overlays rotated boundaries over rotated modulus map and preps for concentric boundary analysis
% % define the x- and y-data for the original line we would like to rotate
% create a matrix of these points, which will be useful in future calculations
v = rot90([yt xt]);
% choose a point which will be the center of rotation
x_center = Xc;
y_center = Yc;
% create a matrix which will be used later in calculations
center = repmat([x_center; y_center], 1, length(xt));
% define angle of rotation
theta = -ang; 
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
% do the rotation...
s = v - center;% shift points in the plane so that the center of rotation is at the origin
vo = R*(v - center) + center; %rotate and shift back to original center
% pick out the vectors of rotated x- and y-data
x_rotated = vo(1,:);
y_rotated = vo(2,:);
% make a plot
modt=mod1(:);
end

