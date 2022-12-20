clear all;
clc;
format long
%% 给定积分面初始流场（单极子问题）
%% give parameter
name1 = 'C:\zzt_ylinder\FW-H-FRE-MOVINGOB-coming-flow\';
outpp_name=[name1,'result_data\dipole\ana_c0=340_M=0.0_numT=4096_motion9-1.dat'];
A=0.02;
rho0=1.0;
c0=340;   %需要注意的是，单极子算例中，声源的波长一定要与积分面的尺寸匹配
u0=0;
f=1;
lamta=f*c0; % characteristic wave length
w=2*pi*f;
dt=0.05;
% set up observer
r(1)=0.0;
r(2)=3400.0;
delta_source=[0.01,0];dxsource=0.02;
M=u0/c0;
beta=(1-M^2)^0.5;
k=w/c0;
%% cal flow
t=0:dt:(4096-1)*dt;
for jt=1:length(t)
%% motion 1 (motion 1-4 are for monopole)
% r(1)=0.0; % 接近声源
% r(2)=34000.0-170*t(jt);
%% motion 2
% r(1)=170*t(jt); % 远离声源
% r(2)=6800.0;
%% motion 3
% r(1)=0; % 周期振荡(垂直)
% r(2)=1360+510*sin(2*pi*0.05*t(jt));
%% motion 4
% r(1)=510*sin(2*pi*0.05*t(jt)); % 周期振荡(水平)
% r(2)=1360;
%% motion 5 (motion 5-6 are for dipole)
% r(1)=-6800+170*t(jt); % 远离声源
% r(2)=6800.0;

%% motion 6
% r(1)=6800+170*t(jt); % 远离声源
% r(2)=0.0;

%% motion 7 (motion 5-6 are for dipole)
% r(1)=-6800+340*0.3*t(jt); % 远离声源
% r(2)=6800.0;


%% motion 8 (motion 5-6 are for dipole)
% r(1)=-6800+340*0.7*t(jt); % 远离声源
% r(2)=6800.0;

%% motion 9 (motion 5-6 are for dipole)
r(1)=-6800+340*1.05*t(jt); % 远离声源
r(2)=6800.0;

% % % % r(1)=5.0*sin(2*pi*2*t(jt)); % 周期转圈
% % % % r(2)=90+5.0*cos(2*pi*2*t(jt));
% % % 
% % % %e_dis=t(end)*2.0; 
% % % %r(1)=0; % 多普勒效应，接近
% % % % r(2)=400+0.5*e_dis-2.0*t(jt);
% % % 
% % % %e_dis=t(end)*2.0; 
% % % %r(1)=0; % 多普勒效应，远离
% % % % r(2)=400-0.5*e_dis+2.0*t(jt);
[ fi1,dfi_dt1,dfi_dx1,dfi_dy1 ] = calfi( t(jt),A,beta,w,M,k,r-delta_source );
% y1轴正向的单极子在（0.01,0）处，因此其对应的r为r-delta_source 
[ fi2,dfi_dt2,dfi_dx2,dfi_dy2 ] = calfi( t(jt),-A,beta,w,M,k,r+delta_source );
% y1轴负向的单极子在（-0.01,0）处，因此其对应的r为r+delta_source 
pp(jt)=-rho0*(dfi_dt1+u0*dfi_dx1+dfi_dt2+u0*dfi_dx2)/dxsource;
up(jt)=(dfi_dx1+dfi_dx2)/dxsource;
vp(jt)=(dfi_dy1+dfi_dy2)/dxsource;
rhop(jt)=pp(jt)/c0^2;
end
%physet0=real([pp]');
%dlmwrite(outpp_name,[t;real([pp])]');
% fft_pp=fft(real(pp));
% fft_pp1=fft_pp';

plot(t,real(pp));
dlmwrite(outpp_name,[t;real([pp])]');
%% determine surface position and output the surface
% 我们用矩形的积分面（以某一点为中心）
%首先确定积分面的参数
pcen(1:2)=[0.0,0.0];
len(1:2)=[0.2*lamta,0.2*lamta];
nump(1:2)=[100,100];
dx=len(1)/nump(1);
dy=len(2)/nump(2);
nodep=zeros(2*nump(1)+2*nump(2),2);

% 计算积分面上的点的位置,由于限制了积分面的形状，我们就不给出每点所在网格的信息
for j=1:nump(1)
    nodep(j,1)=pcen(1)-0.5*len(1)+dx/2.0+(j-1)*dx;
    nodep(j,2)=pcen(2)-0.5*len(2);
    nodep(j+nump(1),1)=nodep(j,1);
    nodep(j+nump(1),2)=pcen(2)+0.5*len(2);
end
for j=1+nump(1)*2:nump(2)+nump(1)*2
    nodep(j,2)=pcen(2)-0.5*len(2)+dy/2.0+(j-1-nump(1)*2)*dy;
    nodep(j,1)=pcen(1)-0.5*len(1);
    nodep(j+nump(2),2)=nodep(j,2);
    nodep(j+nump(2),1)=pcen(1)+0.5*len(1);
end
dlmwrite ([name1,'flow/mesh'],nodep);
%% 输出每个点的速度, 压力，密度随时间变化
for j=1:2*nump(1)+2*nump(2)
%for j=1:2
nameout0=sprintf('flow/flow%d',j+100000);
nameout = [name1,nameout0];
[ fi1,dfi_dt1,dfi_dx1,dfi_dy1 ] = calfi( t,A,beta,w,M,k,nodep(j,:)-delta_source);
[ fi2,dfi_dt2,dfi_dx2,dfi_dy2 ] = calfi( t,-A,beta,w,M,k,nodep(j,:)+delta_source);
pp=-rho0*(dfi_dt1+u0*dfi_dx1+dfi_dt2+u0*dfi_dx2)/dxsource;
up=(dfi_dx1+dfi_dx2)/dxsource;
vp=(dfi_dy1+dfi_dy2)/dxsource;
rhop=pp/c0^2;
physet=real([up;vp;pp;rhop]');
dlmwrite(nameout,physet);    
end












