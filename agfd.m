function [GOR,SEC,Jw,Thb,Tcb,mhn,EE,WR]=agfd (Tcbn,Thb1,L,w,dg,mh1)
N=100; %Number of Intervals for the distance;
N1 =1000;
kg=0.03;%thermal conductivity of air vapor mixture;
dA=L/N*w;%differential area of membrane, m^2;
A=4.6543;%constant in Antonine eqn.;
B=1435.264; %constant in Antonine eqn.;
C=-64.848; %constant in Antonine eqn.;
g=9.8; % gravity constant, kg/s^2
rhol=1000; % water density, kg/m^3
rhov=0.804; % air density, kg/m^3
dp=0.25*10^-3;%thickness of coolant plate, m;
kp=15;%thermal conductivity of coolant plate, W/m/K;

mc=mh1;% coolant mass flow rate, kg/s;
P=101325;% atmospheric pressure, Pa;
Mw=0.018;% molecular weight of water, g/mol;
R=8.314;% universal gas constant, J/mol/K;
kw=0.6;%thermal conductivity of water, W/m/K;

Sa1=3.5;% salinity, wt%;
po_s=0.5;%spacer porosity;
ds=0.8*10^-3;%spacer thickness, m;
df=0.4*10^-3;%spacer filament diameter, m;
dh=4*po_s*df*ds/(2*df+4*(1-po_s)*ds);%hydraulic diameter of spacer, m;;%heat transfer coefficient in feed (W/m^2/K)
Cpl=4180; % specific heat capacity of water, J/kg/K
%------------------
Thb=ones(N,1);
Tf=ones(N,1);
Tcb=ones(N,1);
mh=ones(N,1);
delta=ones(N,1);
hv=ones(N,1);
Sa=ones(N,1);
Q=ones(N,1);
%------------------
mh(1)=mh1; % inlet feed mass flow rate, kg/s;
Thb(1)=Thb1; % inlet feed temperature, K;
Tcb(1)=(Thb1+Tcbn)/2; % guess of outlet coolant temperature, K;
Sa(1)=Sa1; % inlet salinity, wt%;
f=ones(N1,1)*Tcbn;
s=ones(N1,1)*Thb1;
z=0;
Qx=0;

 delta(1)=0;
   for o=1:N1
    for j =1:N
        Tf(j)=(Thb(j)+Tcb(j))/2;
        for i=1:1000000
            for ii=1:1000000
                Th_avg=Thb(j)-273.15;% temperature in degree C
                Sh=Sa(j)*10;% salinity in g/kg
                Cpl_h=(4206.8-6.6197*Sh+1.2288*10^-2*Sh^2) ...
                    +(-1.1262+5.4178*10^-2*Sh-2.2719*10^-4*Sh^2)*Th_avg ...
                    +(1.2026*10^-2-5.3566*10^-4*Sh+1.8906*10^-6*Sh^2)*Th_avg^2 ...
                    +(6.8777*10^-7+1.517*10^-6*Sh-4.4268*10^-9*Sh^2)*Th_avg^3;%water specific heat of feed (J/kg/K)
                mu_Sh=1+(1.474*10^-3+1.5*10^-6*Th_avg-3.927*10^-8*Th_avg^2)*Sh ...
                    +(1.0734*10^-5-8.5*10^-8*Th_avg+2.23*10^-10*Th_avg^2)*Sh^2;%kg/m/s
                mu_h=exp(-3.79418+604.129/(139.18+Th_avg))*mu_Sh*10^-3;%water dynamic viscosity of feed (kg/m/s)
                k_h=10^(log10(240+2*10^-4*Sh) ...
                    +0.434*(2.3-(343+3.7*10^-2*Sh)/(Th_avg+273.15))...
                    *(1-(Th_avg+273.15)/(647.3+3*10^-2*Sh))^(1/3))/1000;%water thermal conductivity of feed (W/m/K)
                Pr_h=Cpl_h*mu_h/k_h;%Prandtl number of feed;
                Re_h=mh(j)/ds/w*dh/mu_h;%local reynold's number in feed
                H_h=0.13*(Re_h)^0.64*Pr_h^0.38*k_h/dh;%heat transfer coefficient in feed (W/m^2/K)
                 %=============================================================================%                     
                a=1.895*10^-5*((Thb(j)+Tf(j))/2)^1.072/P*Mw/R/(dg-delta(j)); %diffusion within air gap
                 %=============================================================================% 
                Jx=a*(10^(A-B/(Thb(j)+C)+5)/(1+5.7357*Sa(j)/(1000-10*Sa(j)))-10^(A-B/(Tf(j)+C)+5));
                hvx=1000*(2501.689845+1.806916015*(Thb(j)-273.15)+ ...
                    5.087717*10^-4*(Thb(j)-273.15)^2-1.1221*10^-5*(Thb(j)-273.15)^3);% enthalpy of saturated water vapor (J/kg)
                Qx1=(Jx*hvx*(dg-delta(j))/kg+Thb(j)-Tf(j))/((dg-delta(j))/kg);% heat flux within feed BL and air gap (J/m^2/s)
                if abs(Qx-Qx1)<1e-4
                    break
                else
                    Qx=(Qx+Qx1)/2;
                    ii=ii+1;
                end
            end
            J(j)=Jx;
            hv(j)=hvx;
            Q(j)=Qx;
            Q_H(j)=J(j)*hv(j)*dA;
            %=============================================================================% 
            Tc_avg=Tcb(j)-273.15;
            Cpl_c=(4206.8)+(-1.1262)*Tc_avg+(1.2026*10^-2)*Tc_avg^2 ...
                +(6.8777*10^-7)*Tc_avg^3;% water specific heat of coolant(J/kg/K)
            mu_c=exp(-3.79418+604.129/(139.18+Tc_avg))*10^-3;% water dynamic viscosity of coolant(kg/m/s);
            k_c=10^(log10(240) ...
                +0.434*(2.3-343/(Tc_avg+273.15))*(1-(Tc_avg+273.15)/647.3)^(1/3));%water thermal conductivity of coolant(W/m/K);
            Re_c=mc/ds/w*dh/mu_c;%local reynold's number in coolant;
            Pr_c=Cpl_c*mu_c/k_c;%Prandtl number;
            H_c=0.13*(Re_c)^0.64*Pr_c^0.38*k_c/dh;%heat transfer coefficient in coolant(W/m^2/K);
            if j==1
            Tf1(j)=Tcb(j)+(delta(j)/kw+dp/kp+1/H_c)*(Q(j)-J(j)*Cpl_c*Tf(j));%temperature at film/plate interface(K)
            else
            Tf1(j)=Tcb(j)+(delta(j)/kw+dp/kp+1/H_c)*(Q(j)-J(j)*Cpl_c*Tf(j)-(mh1-mh(j))*Cpl_c*(Tf(j)-Tf(j-1))/dA);
             end
            if abs(Tf1(j)-Tf(j))<1e-2
                break
            else
                Tf(j)=(Tf(j)+Tf1(j))/2;
                i=i+1;
            end
        end
        mh(j+1)=mh(j)-J(j)*dA;
        Sa(j+1)=mh(j)*Sa(j)/mh(j+1);
        Thb(j+1)=(-Q(j)*dA+mh(j)*Cpl_h*Thb(j))/Cpl_h/mh(j+1);
        if j==1
        Tcb(j+1)=Tcb(j)-(Q(j)-J(j)*Cpl_c*Tf(j))*dA/Cpl_c/mc;
        else
        Tcb(j+1)=Tcb(j)-(Q(j)-J(j)*Cpl_c*Tf(j)-(mh1-mh(j))*Cpl_c*(Tf(j)-Tf(j-1))/dA)*dA/Cpl_c/mc;
        end
        mf(j)=mh1-mh(j);
        mf(j+1)=mh1-mh(j+1);
        Tf_avg=(Thb(j)-Q(j)/H_h+Tf(j))/2-273.15;
        mu_f=exp(-3.79418+604.129/(139.18+Tf_avg))*10^-3;
        v=mu_f/rhol;
        delta(j+1)=(J(j)*3*v/w/g/(rhol-rhov)+delta(j)^3)^(1/3);
    end
    z=Tcb(1);
    if Tcb(N)-Tcbn<0
        f(o)=max(Tcb(1),Tcbn);
    else
        s(o)=min(Tcb(1),Thb1);
    end
    if abs(Tcb(N)-Tcbn)<1e-4
        break

    else
        Tcb(1)=(min(s)+max(f))/2;
    end
    o=o+1;
   end
hv=1000*(2501.689845+1.806916015*(Tcbn-273.15)+ ...
                    5.087717*10^-4*(Tcbn-273.15)^2-1.1221*10^-5*(Tcbn-273.15)^3);
Qp=mh(1)-mh(N);
WR=Qp/mh(1);
Jw=Qp/w/L*3600;
Tcb_1=Tcb(1);
SEC=mh1*Cpl_c*(Thb(1)-Tcb(1))/Qp;
mhn=mh(N);
EE=sum(Q_H)/sum(Q*dA);
GOR=hv/SEC;
delta_N=delta(N);
end