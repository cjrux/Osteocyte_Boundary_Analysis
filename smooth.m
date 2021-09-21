function [finalmat]= smooth(boundmat,mod1,start_size,percent_increase)
for i=1:length(boundmat)
test=boundmat{i}
H=figure(4)
J=plot(test(:,1),test(:,2),'r.-','LineWidth',6)
% xlim([0 length(mod1)])
% ylim([0 length(mod1)])

xlim([-length(mod1) 2*length(mod1)])
ylim([-length(mod1) 2*length(mod1)])

%H=figure(5)
empty=zeros(length(mod1));
im=imagesc(empty);
colormap('white')
hold on;
plot(test(:,1),test(:,2),'k.-','LineWidth',10)
set(gca,'xtick',[])
set(gca,'ytick',[])
%set(gca,'xticklabel',[])
%set(gca,'yticklabel',[])
F=getframe
figure
TestIm=imshow(F.cdata);
[Xim, Mapim] = frame2im(F);
%figure
Testimage=imfill(Xim,'holes');
imshow(Testimage);
bc = strel('diamond',1);
Berode=imerode(Testimage,bc);
JJ=im2bw(Berode);
JJWB=imcomplement(JJ);
Imfilled=imfill(JJWB,'holes');

if i==1
windowSize = start_size; %THIS IS WHERE I'M EDITING
else
    windowSize=round(holdvar*percent_increase)
end
holdvar=windowSize;
kernel = ones(windowSize) / windowSize ^ 2;
blurryImage = conv2(single(Imfilled), kernel, 'same'); 
binaryImage = blurryImage > 0.5; % Rethreshold ORIGINALLY 0.5
imshow([Imfilled, binaryImage]);
%%20
EdgeMap2 = edge(binaryImage,'log');
imshow(EdgeMap2);
%Obtaining Relevant Region Props
s2 = regionprops(EdgeMap2,'PixelList','MajorAxisLength','MinorAxisLength','Centroid');
if size(s2,1)==1
    border2=s2(1).PixelList
else
for k=1:size(s2,1)
placeholder{k}=s2(k).PixelList;%Using that to determine which pixel list defines lacuna edge
end
border2 = cat(1,placeholder{:});
end
% %for k=1:size(s2,1)
% border2=s2(1).PixelList;%Using that to determine which pixel list defines lacuna edge
% %border2=border2
testxxx=border2(:,1)*(length(mod1)/length(Testimage));
testyyy=border2(:,2)*(size(mod1,1)/size(Testimage,1));
% figure
% plot(Xr,Yr,'m','Linewidth',2)
% xlim([0 length(mod1)])
% ylim([0 length(mod1)])
% hold on
% scatter(testxxx,testyyy,'.g','Linewidth',4)

finalmat{i}=[testxxx,testyyy];
close all
end