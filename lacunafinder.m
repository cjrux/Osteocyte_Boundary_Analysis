function [border,mod1,s,I,centroid,EdgeMap]= lacunafinder(mod,element_size,b,a)
mod=rot90(mod)/10^9;%rotating to plot correctly
indices = find(abs(mod)<5);%filtering out data less than 5 GPa
mod(indices)=0; %replacing points with NaN
figure(1)
imagesc(mod)
colormap(hot)
 hcb = colorbar;
 title(hcb,'GPa')
    %a=30; %colorbar upper bound usually 30
   %b=0; %colorbar lower bound
    caxis([b a])%these values can be changed depending on modulus range
    pbaspect([1 1 1]);
%---Guassian Filter---
mod1 = imgaussfilt(mod, 1);
mod1(isnan(mod1))=0;
new_map = mod1(:,:) > median(median(mod1)) - 3*(std(std(mod1)));%may need changes #3* is standard
%---Spot Correction---
spot_correction = bwareaopen(imcomplement(new_map), 1000); %originally 1000
%---Sphere Fit Erosion---
%ab = strel('sphere',5);%sphere erosion with dictated size -> originally 2
ab = strel('diamond',element_size);%originally 5
erodedBWs = imerode(spot_correction,ab);
EdgeMap = edge(erodedBWs,'log');
%Obtaining Relevant Region Props
s = regionprops(EdgeMap,'PixelList','MajorAxisLength','MinorAxisLength','Centroid','Orientation','Area');%obtaining pixel list of all boundary objects
A=cell2mat(struct2cell(regionprops(EdgeMap,'MajorAxisLength')));%obtaining all major axis lengths for objects
I=find(A == max(A(:)));%finding location of the max Major Axis Length
border=s(I).PixelList;%Using that to determine which pixel list defines lacuna edge
centroid=s(I).Centroid;
end