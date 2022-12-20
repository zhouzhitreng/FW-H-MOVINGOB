clear all;
clc;
format long
%% �����������ʼ���������������⣩
%% give parameter
outpp_name='../../result_data/monopole/ana_c0=340_M=0.0_numT=2048_motion3.dat';
A=0.02;
rho0=1.0;
c0=340;   %��Ҫע����ǣ������������У���Դ�Ĳ���һ��Ҫ�������ĳߴ�ƥ��
u0=0.0;
f=1;
lamta=f*c0; % characteristic wave length
w=2*pi*f;
dt=0.05;
% set up observer
r(1)=0.0;
r(2)=3400.0;
M=u0/c0;
beta=(1-M^2)^0.5;
k=w/c0;
%% cal flow
t=0:dt:(2048-1)*dt;
for jt=1:length(t)
%% motion 1
% r(1)=0.0; % �ӽ���Դ
% r(2)=34000.0-170*t(jt);
%% motion 2
% r(1)=170*t(jt); % Զ����Դ
% r(2)=6800.0;
%% motion 3
r(1)=0; % ������(��ֱ)
r(2)=1360+510*sin(2*pi*0.05*t(jt));
%% motion 4
% r(1)=510*sin(2*pi*0.05*t(jt)); % ������(ˮƽ)
% r(2)=1360;



% % % % r(1)=5.0*sin(2*pi*2*t(jt)); % ����תȦ
% % % % r(2)=90+5.0*cos(2*pi*2*t(jt));
% % % 
% % % %e_dis=t(end)*2.0; 
% % % %r(1)=0; % ������ЧӦ���ӽ�
% % % % r(2)=400+0.5*e_dis-2.0*t(jt);
% % % 
% % % %e_dis=t(end)*2.0; 
% % % %r(1)=0; % ������ЧӦ��Զ��
% % % % r(2)=400-0.5*e_dis+2.0*t(jt);
[ fi,dfi_dt,dfi_dx,dfi_dy ] = calfi( t(jt),A,beta,w,M,k,r );
pp(jt)=-rho0*(dfi_dt+u0*dfi_dx);
up(jt)=dfi_dx;
vp(jt)=dfi_dy;
rhop(jt)=pp(jt)/c0^2;
end
%physet0=real([pp]');
%dlmwrite(outpp_name,[t;real([pp])]');
% fft_pp=fft(real(pp));
% fft_pp1=fft_pp';

plot(t,real(pp));
dlmwrite(outpp_name,[t;real([pp])]');
%% determine surface position and output the surface
% �����þ��εĻ����棨��ĳһ��Ϊ���ģ�
%����ȷ��������Ĳ���
pcen(1:2)=[0.0,0.0];
len(1:2)=[0.2*lamta,0.2*lamta];
nump(1:2)=[100,100];
dx=len(1)/nump(1);
dy=len(2)/nump(2);
nodep=zeros(2*nump(1)+2*nump(2),2);

% ����������ϵĵ��λ��,���������˻��������״�����ǾͲ�����ÿ�������������Ϣ
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
dlmwrite ('../../flow/mesh',nodep);
%% ���ÿ������ٶ�, ѹ�����ܶ���ʱ��仯
for j=1:2*nump(1)+2*nump(2)
%for j=1:2
nameout=sprintf('../../flow/flow%d',j+100000);
 [ fi,dfi_dt,dfi_dx,dfi_dy ] = calfi( t,A,beta,w,M,k,nodep(j,:));
pp=-rho0*(dfi_dt+u0*dfi_dx);
up=dfi_dx;
vp=dfi_dy;
rhop=pp/c0^2;
physet=real([up;vp;pp;rhop]');
dlmwrite(nameout,physet);    
end












