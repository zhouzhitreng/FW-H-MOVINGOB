function [ green,dg_dkesi,dg_dita ] = calgreen( beta,M,k,r )
% use numercial method for difference
format long
r0=(r(1)^2+beta^2*r(2)^2)^0.5;
kk1=abs(k);
hankel=besselh(0,2,r0*kk1/beta^2);
dhankel=(besselh(0,2,r0*kk1/beta^2+0.00001)-besselh(0,2,r0*kk1/beta^2-0.00001))/0.00002;
if (k<0)
    hankel=-hankel;
    dhankel=dhankel;
end
green=(i/(4.0*beta))*exp(i*M*k*r(1)/beta^2)*hankel;
dg_dkesi1=-i*M*k*green/beta^2;
dg_dkesi2=(i/(4.0*beta))*exp(i*M*k*r(1)/beta^2)*dhankel*k*(-r(1))/(beta^2*r0);
dg_dkesi=dg_dkesi1+dg_dkesi2;
dg_dita=(i/(4.0*beta))*exp(i*M*k*r(1)/beta^2)*dhankel*k*(-r(2))/r0;

end

