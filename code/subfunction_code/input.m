clc;
clear all;
% number of points
nx=200;
ny=200;
ncore=100;
%%%%%%%%%%%%%%%%%%
%position of start point
xst=-8/2;
yst=-8/2;
%%%%%%%%%%%%%%%%%%
% number of cores
% number of points in each core
%nco=200;
npp=nx*ny/ncore;

for i=1:nx
    for j=1:ny
    a((i-1)*nx+j,1)=xst+(j-1)*2*abs(xst)/199;
    a((i-1)*nx+j,2)=yst+(i-1)*2*abs(xst)/199;
    a((i-1)*nx+j,3)=0;
    end   
end

for i=0:ncore-1
%dirs=dir(['input',i])
aaa=sprintf('input%d',i);
%dlmwrite(aaa,a(i*npp+1:i*npp+npp,:));
end

%% prepare points for calculating directivity
dtheta=2.0
nspl=floor(360.0/dtheta);
for i=1:nspl
   xx(i)=10.0*cos((i-1)*dtheta*3.141592653/180);
   yy(i)=10.0*sin((i-1)*dtheta*3.141592653/180);
   zz(i)=0.0
end
aspl=[xx;yy;zz]'
dlmwrite('directivity',aspl);








