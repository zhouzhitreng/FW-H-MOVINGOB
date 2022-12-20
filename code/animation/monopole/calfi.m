function [ fi,dfi_dt,dfi_dx,dfi_dy ] = calfi( t,A,beta,w,M,k,r )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
r0=(r(1)^2+beta^2*r(2)^2)^0.5;
%kk1=abs(k);
hankel=besselh(0,2,r0*k/beta^2);
dhankel=(besselh(0,2,r0*k/beta^2+0.01)-besselh(0,2,r0*k/beta^2-0.01))/0.02;
if (k<0)
   dhankel=-dhankel;
end
% hankel=real(besselj(0,r0*k/beta^2))-i*real(bessely(0,r0*k/beta^2));
% dhankel1=(real(besselj(0,r0*k/beta^2+0.00001))-real(besselj(0,r0*k/beta^2-0.00001)))/0.00002;
% dhankel2=(real(bessely(0,r0*k/beta^2+0.00001))-real(bessely(0,r0*k/beta^2-0.00001)))/0.00002;
% dhankel=dhankel1-i*dhankel2;
% 差分法求导数这一句我不知道有没有问题，应该没有
fi=A*(i/(4.0*beta))*exp(i*(w*t+M*k*r(1)/beta^2))*hankel;
dfi_dt=i*w*fi;
dfi_dx1=i*M*k*fi/beta^2;
dfi_dx2=A*(i/(4.0*beta))*exp(i*(w*t+M*k*r(1)/beta^2))*dhankel*k*(r(1))/(beta^2*r0);
dfi_dx=dfi_dx1+dfi_dx2;
dfi_dy=A*(i/(4.0*beta))*exp(i*(w*t+M*k*r(1)/beta^2))*dhankel*k*(r(2))/r0;

end

