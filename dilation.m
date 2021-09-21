function boundmat = dilation(border, major, minor,c,num_dilations,mod,map_dim,step_size)
for k=1:num_dilations

for i=1:length(border(:,1))
    xvec(i)=border(i,1)-c(1,1)+2*(major/minor)*(border(i,1)-c(1,1)/(max(border(:,1)))-c(1,1));
    yvec(i)=border(i,2)-c(1,2);
    normx(i)=xvec(i)/sqrt(xvec(i)^2+yvec(i)^2);
    normy(i)=yvec(i)/sqrt(xvec(i)^2+yvec(i)^2);
    
end
unit_V(:,1)=normx;
unit_V(:,2)=normy;


Scaled_V=k*unit_V*step_size*length(mod)/map_dim;%Here's where to change scale of dilation initially 320/12


new_bound=Scaled_V+border;
boundmat{k}=new_bound;
end
