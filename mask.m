function [modx,inmat]= mask(mod,element_size,lowerbound,upperbound)
%% Actual Lacunar Border -- points used for actual appropriately eroded lacuna edge
%best to capture whole osteocyte, typically 4 -- determines intensity of erosion
[lacunaredge,mod2,s2,I2,centroid,EdgeMap]=lacunafinder(mod,element_size,lowerbound,upperbound); %image processing function
hold on
scatter(lacunaredge(:,1),lacunaredge(:,2),'c.') %can check border overlayed on mod map here
[x_lacunar_edge,y_lacunar_edge]= CCW(centroid, lacunaredge);
%% 
n=length(mod); %determining size of matrix 
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

%%
hold off
in = inpolygon(xt,yt,x_lacunar_edge,y_lacunar_edge);
in=~in;
in_num=double(in);
inmat=reshape(in_num,length(mod),length(mod));
%imagesc(inmat)
%%
modx=(rot90(mod).*(inmat));
%%
%clims = [0 1];
%imagesc(modx, clims)