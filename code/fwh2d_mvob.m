clear all;
clc;
format long
%% 我只接受水平来流
%认为远场没有直流的信息
% input parameter
% physical parameter
%outpp_name='D:\zhouzhiteng\research\ACOUSTICS\code\2dFWH_copy_1\0_result\dt_test\1_monopole_num_U=0.0_dt=0010.dat';
example='dipole';
adpath(example);
outpp_name=['result_data/',example,'/num_c0=340_M=0.0_numT=4096_motion5.dat'];
T=1.0;
rho0=1.0d0;
p0=0.0;
c0=340;
lamta=340;  % characteristic wave length
u0=0;
M=abs(u0)/c0;
beta=(1-M^2)^0.5;
%parameter about output time
dt=0.05;
Fs=1/dt;
t=0:dt:4095*dt;
Len=length(t);
nt1_final=Len;nt1_old=nt1_final;
Len_half=(Len-1)/2;
pp=zeros(1,Len);
pp1=pp;
pp2=pp;
pp3=pp;
pq=zeros(1,Len); %pq is for test
%geometry information
%矩形的四条边按照下上左右的顺序安排
npoint=90;
%prms=zeros(1,90);
%theta=zeros(1,90);
%for inp=1:npoint
pp=zeros(1,Len);
pp1=pp;
pp2=pp;
pp3=pp;
 %   theta(inp)=(inp-1)*360.0/(npoint*1.0);
%ob_pos=[c0*T*cos(2*pi*theta(inp)/360.0),c0*T*sin(2*pi*theta(inp)/360.0)];
%ob_pos=[75*cos(2*pi*theta(inp)/360.0),75*sin(2*pi*theta(inp)/360.0)];
%ob_pos=[0,90];
len(1:2)=[0.2*lamta,0.2*lamta];
nump(1:2)=[100,100];
dx=len(1)/nump(1);
dy=len(2)/nump(2);
mesh=dlmread('flow/mesh'); %mesh is the position of nodes on the surface
%%  set up moving observer
%% motion 1 (motion 1-4 are for monopole)
% ob_mvpos(:,1)=0.0*t;
% ob_mvpos(:,2)=34000-170*t;

%% motion 2
% ob_mvpos(:,1)=170*t; % 接近声源
% ob_mvpos(:,2)=6800+0.0*t;

%% motion3
% ob_mvpos(:,1)=0*t;% 周期振荡(垂直)
% ob_mvpos(:,2)=1360+510*sin(2*pi*0.05*t);

%% motion4
% ob_mvpos(:,1)=510*sin(2*pi*0.05*t);% 周期振荡(水平)
% ob_mvpos(:,2)=1360;

%% motion 5 (motion 5-6 are for dipole)
% we can inply that the moment crossing through the ... is 40s
ob_mvpos(:,1)=-6800+170.0*t; % 接近声源
ob_mvpos(:,2)=0.0*t+6800;
 

%% motion 6
% ob_mvpos(:,1)=6800+170.0*t; % 接近声源
% ob_mvpos(:,2)=0.0*t;

%% motion 7 (motion 5-6 are for dipole)
% ob_mvpos(:,1)=-6800+340*0.3*t; % 接近声源
% ob_mvpos(:,2)=0.0*t+6800;

%% motion 8 (motion 5-6 are for dipole)
% ob_mvpos(:,1)=-6800+340*0.7*t; % 接近声源
% ob_mvpos(:,2)=0.0*t+6800;

%% motion 9 (motion 5-6 are for dipole)
% ob_mvpos(:,1)=-6800+340*1.05*t; 
% ob_mvpos(:,2)=0.0*t+6800;



% % % % ob_mvpos(:,1)=5.0*sin(2*pi*2.0*t);% 周期转圈
% % % % ob_mvpos(:,2)=90+5.0*cos(2*pi*2.0*t);
% % % 
% % % %e_dis=t(end)*2.0; 
% % % %ob_mvpos(:,1)=0; % 多普勒效应，接近
% % % %ob_mvpos(:,2)=400+0.5*e_dis-2.0*t;
% % % 
% % % %e_dis=t(end)*2.0; 
% % % %ob_mvpos(:,1)=0; % 多普勒效应，远离
% % % %ob_mvpos(:,2)=400-0.5*e_dis+2.0*t;
%%  cal source in time domain
n=[0,-1];
dl=dx;
for j=1:2*nump(1)+2*nump(2)
  
%for j=2*nump(1)+nump(2)+1:2*nump(1)+nump(2)+1
%for j=1:1
    r=((ob_mvpos(:,1)-mesh(j,1)).^2+beta^2*(ob_mvpos(:,2)-mesh(j,2)).^2).^0.5;
    t1_new=t'-(r+M*mesh(j,1)-M*ob_mvpos(:,1))/c0/beta^2;
    t1_new=t1_new';
    if (j==(nump(1)+1))
        n=[0,1];
    end
    if (j==(2*nump(1)+1))
        n=[-1,0];
        dl=dy;
    end
    if (j==(2*nump(1)+nump(2)+1))
        n=[1,0];
    end
namein=sprintf('flow/flow%d',j+100000);
flow_var1=dlmread(namein);
flow_var1(:,1:2)=flow_var1(:,1:2);
flow_var1(:,1)=flow_var1(:,1);%-u0;  %!zzt
flow_var1(:,3)=flow_var1(:,3);

flow_var(:,1:3)=flow_var1(:,1:3);
flow_var(:,4)=0.0;
%for testm=1:10
%    flow_var(1+(testm-1)*99:1+(testm-1)*99+98,1:3)=flow_var1(33:131,:);
%end
%flow_var is the set of variables(up,vp,pp,rhop);
p=p0+flow_var(:,3);
flow_var(:,4)=flow_var(:,3)/c0^2;
rho=rho0+flow_var(:,4);
tempnu0=n(1)*(flow_var(:,1))+n(2)*(flow_var(:,2));
tempnu=tempnu0+n(1)*u0;
%F1=n(1)*p+diag(rho*(diag((flow_var(:,1)-u0)*tempnu'))')+rho0*n(1)*u0^2;
F1=n(1)*p+((rho.*(flow_var(:,1)-u0)).*tempnu)+rho0*n(1)*u0^2;
temp0=flow_var(:,1)-u0;
temp1=temp0.*tempnu;
temp2=rho*temp1';

F2=n(2)*p+(rho.*(((flow_var(:,2)-0.0).*tempnu)));
Q=rho.*tempnu-rho0*n(1)*u0;
%% conduct fft on source important!!!!!
FFT_F1=fft(F1);
FFT_F2=fft(F2);
FFT_Q=fft(Q);
fre=Fs*(0:(Len/2))/Len;
for frei=2:length(fre)
    w1=2*pi*fre(frei);  
    k1=w1/c0;   
    %[ green,dg_dkesi,dg_dita ] = calgreen( beta,M,k1,r );
    coff_F=(1.0/4.0)*(2.0/pi/k1)^0.5*w1/c0/beta^2;
    coff_Q=(-i/4.0)*(2.0/pi/k1)^0.5;
    pp1(1,frei)=coff_F*FFT_F1(frei)*exp(i*pi/4.0);
    pp2(1,frei)=coff_F*FFT_F2(frei)*exp(i*pi/4.0);
    pp3(1,frei)=coff_Q*i*w1*FFT_Q(frei)*exp(i*pi/4.0);
    
    frei1=Len-frei+2;
    pp1(1,frei1)=conj(pp1(1,frei));
    pp2(1,frei1)=conj(pp2(1,frei));
    pp3(1,frei1)=conj(pp3(1,frei));
    %pp(1,frei1)=conj(pp(1,frei));
    pp1p=pp1';
    pp2p=pp2';
    pp3p=pp3';
end
% ifft is needed for each point three times(F1,F2,Q)
ifft_F1=real(ifft(pp1));ifft_F2=real(ifft(pp2));
ifft_Q=real(ifft(pp3));
% intepolate or not????
[t1_new_arr,t1_old_arr,nt1]=pre_intep(t1_new,t);
[ifft_F1_part,ifft_F2_part,ifft_Q_part]=intep_part(ifft_F1,ifft_F2,ifft_Q,t1_new_arr,t1_old_arr,nt1);
% terms else
ifft_F1=dl*ifft_F1_part.*(1.0./r(end-nt1+1:end)').^0.5.*...
                                    (-M+(ob_mvpos(end-nt1+1:end,1)'-mesh(j,1))./r(end-nt1+1:end)');
ifft_F2=dl*ifft_F2_part.*(1.0./r(end-nt1+1:end)').^0.5.*...
                                    (ob_mvpos(end-nt1+1:end,2)'-mesh(j,2))*beta^2./r(end-nt1+1:end)';
ifft_Q=dl*ifft_Q_part.*(1.0./r(end-nt1+1:end)').^0.5;
% add the present point to all
nt1_final=min([nt1,nt1_old]);
pp(end-nt1_final+1:end)=pp(end-nt1_final+1:end)+...
                                ifft_Q(end-nt1_final+1:end)+...
                                ifft_F1(end-nt1_final+1:end)+...
                                ifft_F2(end-nt1_final+1:end);
nt1_old=nt1_final;
end
plot(t,real(pp));
dlmwrite(outpp_name,[t;real([pp])]');
%% conduct ifft
%pp_t=ifft(pp);
%ppp=pp';
%hold on;
%plot(t,pp_t);
% %dlmwrite(outpp_name,[t;pp_t]');
% prms(1,inp)=0.7*max(pp_t(200:300));
% end
% dlmwrite('prms.dat',real([theta;prms]'));
