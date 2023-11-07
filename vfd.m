function [GOR,SEC,Jw,Thb,WR]=vfd(Thb1,L,w,mh1,Pv)
N=100; %Number of Intervals for the distance;
N1 =1000;

dg=100*10^-6;%membrane thickness, m;

kg=0.03;%thermal conductivity of air vapor mixture;

dA=L/N*w;%differential area of membrane, m^2;
A=4.6543;%constant in Antonine eqn.;
B=1435.264; %constant in Antonine eqn.;
C=-64.848; %constant in Antonine eqn.;
g=9.8; % gravity constant, kg/s^2
rhol=1000; % water density, kg/m^3
rhov=0.804; % air density, kg/m^3


P=101325;% atmospheric pressure, Pa;
Mw=0.018;% molecular weight of water, g/mol;
R=8.314;% universal gas constant, J/mol/K;

kw=0.6;%thermal conductivity of water, W/m/K;

Sa1=3.5;% salinity, wt%;
po_s=0.5;%spacer porosity;
ds=0.8*10^-3;%spacer thickness, m;
df=0.4*10^-3;%spacer filament diameter, m;
dh=4*po_s*df*ds/(2*df+4*(1-po_s)*ds);%hydraulic diameter of spacer, m;
Cpl=4180; % specific heat capacity of water, J/kg/K

%------------------
Thb=ones(N,1);
Thm=ones(N,1);
mh=ones(N,1);
hv=ones(N,1);
Sa=ones(N,1);
Q=ones(N,1);
H_h=ones(N,1);
RA=ones(N,1);
%------------------
mh(1)=mh1;
Thb(1)=Thb1;
Sa(1)=Sa1;
z=0;
Qx=0;
To=298;

    for j =1:N
            for ii=1:1000000
                Th_avg=Thb(j)-273.15;
                Sh=Sa(j)*10;% salinity in g/kg
                Cpl_h=(4206.8-6.6197*Sh+1.2288*10^-2*Sh^2) ...
                    +(-1.1262+5.4178*10^-2*Sh-2.2719*10^-4*Sh^2)*Th_avg ...
                    +(1.2026*10^-2-5.3566*10^-4*Sh+1.8906*10^-6*Sh^2)*Th_avg^2 ...
                    +(6.8777*10^-7+1.517*10^-6*Sh-4.4268*10^-9*Sh^2)*Th_avg^3;%water specific heat of feed(J/kg/K)
                mu_Sh=1+(1.474*10^-3+1.5*10^-6*Th_avg-3.927*10^-8*Th_avg^2)*Sh ...
                    +(1.0734*10^-5-8.5*10^-8*Th_avg+2.23*10^-10*Th_avg^2)*Sh^2;%kg/m/s
                mu_h=exp(-3.79418+604.129/(139.18+Th_avg))*mu_Sh*10^-3;%water dynamic viscosity of feed(kg/m/s)
                k_h=10^(log10(240+2*10^-4*Sh) ...
                    +0.434*(2.3-(343+3.7*10^-2*Sh)/(Th_avg+273.15))...
                    *(1-(Th_avg+273.15)/(647.3+3*10^-2*Sh))^(1/3));%water thermal conductivity of feed(W/m/K)
                Pr_h=Cpl_h*mu_h/k_h;%Prandtl number of feed;
                Re_h=mh(j)/ds/w*dh/mu_h;%local reynold's number in feed
                H_h=0.13*(Re_h)^0.64*Pr_h^0.38*k_h/dh;%heat transfer coefficient in feed (W/m^2/K)
                
                 %=============================================================================%
                a=1.895*10^-5*(Thb(j)-Qx/H_h)^1.072/P*Mw/R/dg; %diffusion within gap                  
                 %=============================================================================% 
                Jx=a*(10^(A-B/(Thb(j)-Qx/H_h+C)+5)/(1+5.7357*Sa(j)/(1000-10*Sa(j)))-Pv);
                hvx=1000*(2501.689845+1.806916015*(Thb(j)-Qx/H_h-273.15)+ ...
                    5.087717*10^-4*(Thb(j)-Qx/H_h-273.15)^2-1.1221*10^-5*(Thb(j)-Qx/H_h-273.15)^3);% enthalpy of saturated water vapor (J/kg)
                Qx1=Jx*hvx;% heat flux (J/m^2/s)
                if abs(Qx-Qx1)<1e-5
                    break
                else
                    Qx=(Qx+Qx1)/2;
                    ii=ii+1;
                end
            end
        J(j)=Jx;
        hv(j)=hvx;
        Q(j)=Qx;
        RA(j)=1/a;
        mh(j+1)=mh(j)-J(j)*dA;
        Sa(j+1)=mh(j)*Sa(j)/mh(j+1);
        Thb(j+1)=(-Q(j)*dA+mh(j)*Cpl_h*(Thb(j)-273.15))/Cpl_h/mh(j+1)+273.15;
        end

hv=1000*(2501.689845+1.806916015*(Thb(N)-273.15)+ ...
                    5.087717*10^-4*(Thb(N)-273.15)^2-1.1221*10^-5*(Thb(N)-273.15)^3);
Qp=mh(1)-mh(N);
Jw=Qp/w/L*3600;
WR=Qp/mh(1)*100;
SEC=mh(1)*Cpl_h*(Thb(1)-To)/Qp;
mhn=mh(N);
GOR=hv/SEC;
end


