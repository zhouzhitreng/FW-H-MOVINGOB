clear all;
clc;
format long
%% 给定积分面初始流场（单极子问题）
%% give parameter
outpp_name='tec_c0340_u0_170_t1.0.dat';
A=0.02;
rho0=1.0;
c0=340;   %需要注意的是，单极子算例中，声源的波长一定要与积分面的尺寸匹配
u0=170.0;
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
%% 绘制云图
nump(1)=680;nump(2)=680;
dx=10;dy=10;
len(1)=nump(1)*dx;
len(2)=nump(2)*dy;
t_plot=1.0;
for j2=1:nump(2)
    for j1=1:nump(1)
        j=j1+(j2-1)*nump(1);
        node_plot(j,1)=dx/2.0-(len(1)+dx)/2.0+(j1-1)*dx;
        node_plot(j,2)=dy/2.0-(len(2)+dy)/2.0+(j2-1)*dy;
        [ fi,dfi_dt,dfi_dx,dfi_dy ] = calfi( t_plot,A,beta,w,M,k,node_plot(j,:));
        pp_plot(j)=-rho0*(dfi_dt+u0*dfi_dx);
        up_plot(j)=dfi_dx;
        vp_plot(j)=dfi_dy;
    end
end

fid=fopen(outpp_name,'w');
    fid=fopen(outpp_name,'w');
    fprintf(fid,'TITLE  ="AMR2D"\nVARIABLES = "x"\n"y"\n"p"\n"u"\n"v"\nZONE T="ZONE T=BLKPROC000001"\n STRANDID=0, SOLUTIONTIME=0\n');
    fprintf(fid,'I=%d J=%d ZONETYPE=Ordered\n',nump(1),nump(2));
    fprintf(fid,'DATAPACKING=POINT\n');
    fprintf(fid,'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)\n');
    for l=1:nump(1)*nump(2);
        if (node_plot(l,1)==0&&node_plot(l,2)==0) 
            pp_plot(l)=0;up_plot(l)=0;vp_plot(j)=0;
        end
       fprintf(fid,'%d %d %d %d %d\n',node_plot(l,1),node_plot(l,2),pp_plot(l),up_plot(l),vp_plot(j));
    end
fclose(fid);










