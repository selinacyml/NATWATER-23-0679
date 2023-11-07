function [GOR,SEC,Jw,Thb,Ta,WR]=sgfd(Ta1,Thbn,L,w,mhn,ma1)
N=100; %Number of Intervals for the distance;
N1 =10000;
dg=100*10^-6;%membrane thickness, m;
t=2;%tortuosity;
po=0.75;% porosity;
kg=0.028;%thermal conductivity of air vapor mixture;
dA=L/N*w;%differential area of membrane, m^2;
A=4.6543;% constant in Antonine eqn.;
B=1435.264; % constant in Antonine eqn.;
C=-64.848; % constant in Antonine eqn.;
g=9.8; % gravity accelaration m/s^2
rhol=100; % water density kg/m^3
rhov=0.804; % air density kg/m^3

P=101325;% atmospheric pressure, Pa;
Ptot=ones(N,1)*P; % total pressure in the air channel, Pa (initial values)
Mw=0.018;% molecular weight of water, kg/mol;
Ma=0.029;% molecular weight of air, kg/mol;
R=8.314;% universal gas constant, J/mol/K; Pa*m3/mol/K
k=0.18*(1-po)+0.028*po; % thermal conductivity of membrane (W/m/K);
kw=0.6;%thermal conductivity of water, W/m/K;

Sa1=0; %salinity, wt%;
po_s=0.5;% spacer porosity;
%ds=0.002;%spacer thickness, m
dh=0.002;%hydraulic diameter of spacer, m;
d_ch=0.1; % air channel thickness, m
dha=min(2*d_ch*w/(w+d_ch),d_ch); % hydraulic diameter of channel, m 
Cpl=4180; % specific heat capacity of water, J/kg/K

ds=0.8*10^-3;%spacer thickness, m;
df=0.4*10^-3;%spacer filament diameter, m;
dhf=4*po_s*df*ds/(2*df+4*(1-po_s)*ds);%hydraulic diameter of spacer, m;


%------------------
mh=mhn*ones(N,1); % initial values of mass flow rate of the feed in the feed channel, kg/s
Sa=zeros(N,1); % initial values of salinity
m_a=ones(N,1)*ma1;
mh0=0;

mh_s=ones(N1,1)*mhn;
mh_b=ones(N1,1)*(mhn+mhn);

%------------------
m_a(1)=ma1; % mass flow rate of feed at the inlet, kg/s
Ta(1)=Ta1; % feed temperature of the feed at the inlet, K
%(Thbn+Ta1)/2; % guess of initial value of the coolant at the outlet, K
Sa(1)=Sa1; % feed salinity at the inlet
q_a=ma1/Ma*R*Ta1/P; % volume flow rate of air at the inlet, m^3/s
%P*q_a/R/Tan*Ma; % mass flow rate of air at the inlet, kg/s
Cvw=1460; % heat capacity of vapor, J/kg/K
Cva=1000; % heat capacity of air, J/kg/K
m_w=zeros(N,1); % mass flow rate of vapor in the channel, kg/s
rh=0;
m_w0=10^(A-B/(Ta1+C)+5)/P*ma1/29*18*rh; % mass flow rate of vapor in the channel at the inlet, kg/s
m_w(1)=m_w0;
rho=1000;
rhoa=1.3;
To=298;

v=ma1/rhoa/d_ch/w;

% Q1=1000000;
% Q2=0;
% Qx=(Q1+Q2)/2;

Thb1=Thbn;
Thb2=Ta1;
Thb(1)=(Thb1+Thb2)/2;
%------------------------------------------------------
for oo=1:N1
    mh1=mhn;
    mh2=0;
    mh(1)=(mh1+mh2)/2;
    
    for o=1:N1
        for j =1:N-1
            Q1=100000000;
            Q2=0;
            Qx=(Q1+Q2)/2;
            for ii=1:1000000
                Th_avg=Thb(j)-273.15;
                Sh=Sa(j)*10;
                Cpl_h=(4206.8-6.6197*Sh+1.2288*10^-2*Sh^2) ...
                    +(-1.1262+5.4178*10^-2*Sh-2.2719*10^-4*Sh^2)*Th_avg ...
                    +(1.2026*10^-2-5.3566*10^-4*Sh+1.8906*10^-6*Sh^2)*Th_avg^2 ...
                    +(6.8777*10^-7+1.517*10^-6*Sh-4.4268*10^-9*Sh^2)*Th_avg^3;
                mu_Sh=1+(1.474*10^-3+1.5*10^-6*Th_avg-3.927*10^-8*Th_avg^2)*Sh ...
                    +(1.0734*10^-5-8.5*10^-8*Th_avg+2.23*10^-10*Th_avg^2)*Sh^2;
                mu_h=exp(-3.79418+604.129/(139.18+Th_avg))*mu_Sh*10^-3;    
                k_h=10^(log10(240+2*10^-4*Sh) ...
                    +0.434*(2.3-(343+3.7*10^-2*Sh)/(Th_avg+273.15))...
                    *(1-(Th_avg+273.15)/(647.3+3*10^-2*Sh))^(1/3));%water thermal conductivity of feed(W/m/K)
                Pr_h=Cpl_h*mu_h/k_h;
                Re_h=mh(j)/ds/w*dhf/mu_h;%local reynold's number in feed  
                H_h=0.13*(Re_h)^0.64*Pr_h^0.38*k_h/dhf;%heat transfer coefficient in feed (W/m^2/K)
                Nu_h=0.664*Re_h^0.5*Pr_h^0.33;% convective heat transfer coefficient of air, W/m2/K
               %=============================================================================%
               a=1.895*10^-5*((Thb(j)-Qx/H_h+Ta(j))/2)^1.072/P*Mw/R/dg; %diffusion within gap
              %=============================================================================%  
              Jx=a*(10^(A-B/(Thb(j)+Qx/H_h+C)+5)-Ptot(j)*(m_w(j)/Mw)/(m_w(j)/Mw+ma1/Ma));   
              hvx=1000*(2501.689845+1.806916015*(Thb(j)+Qx/H_h-273.15)+ ...
              5.087717*10^-4*(Thb(j)+Qx/H_h-273.15)^2-1.1221*10^-5*(Thb(j)+Qx/H_h-273.15)^3);% enthalpy of saturated water vapor, J/kg
              %=============================================================================%        
              Ta_avg=Ta(j);
              Cpa=1000*(1.03409-0.284887*10^-3*Ta_avg+0.7816818*10^-6*Ta_avg^2-0.4970786*10^-9*Ta_avg^3+0.1077024*10^-12*Ta_avg^4); % heat capacity of air, J/kg/K
              mu_a=10^-6*(-0.98601+9.080123*10^-2*Ta_avg-1.17635575*10^-4*Ta_avg^2+1.2349703*10^-1*Ta_avg^3-5.7971299*10^-11*Ta_avg^4); % viscosity of air, Ns/m2
              k_a=-2.276501*10^-6+1.2598485*10^-7*Ta_avg-1.4815235*10^-10*Ta_avg^2+1.73550646*10^-13*Ta_avg^3 ...
              -1.066657*10^-16*Ta_avg^4+2.47663035*10^-20*Ta_avg^5; % air thermal conductivity, W/m/K
              Pr_a=Cpa*mu_a/k_a;
              Re_a=ma1/Ma/d_ch/w*dha/mu_a;%local reynold's number in air
              Nu_a=0.664*Re_a^0.5*Pr_a^0.33;
              H_a=Nu_a*k_a/dha;% convective heat transfer coefficient of air, W/m2/K
              %=============================================================================%
        
              Qx1=((Thb(j)-Ta(j))*(1+Jx*Cpv)+Jx*hvx*dg/kg+Jx*hvx/H_a)/(1/H_h+dg/kg+1/H_a+Jx*Cpv/H_a/H_h);
  
               if abs((Qx1-Qx)/Qx) < 0.01 
                    break
                else
                     Qx=(Qx+Qx1)/2;
                end
               
            ii=ii+1;
 
            end
            J(j)=Jx;
            Q(j)=Qx;
            RMa(j)=1/a;
            m_w(j+1)=m_w(j)+J(j)*dA;
            Cpv=1000*(1.86910989-2.578421578*10^-4*Tma(j)+1.941058941*10^-5*Tma(j)^2);
            Ta(j+1)=((m_a(j)*Cva+m_w(j)*Cpv)*(Ta(j)-273.15)+(Qx-Jx*hvx+Jx*Cpv*(Thb(j)-Ta(j)-Qx/H_h))*dA)/(m_a(j)*Cva+m_w(j+1)*Cpv)+273.15;
            mh(j+1)=mh(j)+J(j)*dA;
            Thb(j+1)=(Qx*dA+mh(j)*Cpl_h*(Thb(j)-273.15))/Cpl_h/mh(j+1)+273.15;
        end
        
      
        if abs(mh(N)-mhn)< mhn/100000
            break
        elseif mh(N)-mhn<0
                  mh2=max(mh(1),mh2);
              else
                  mh1=min(mh(1),mh1);
        end
                   mh(1)=(mh1+mh2)/2
    
        o=o+1;
   end
   


        if abs(Thb(N)-Thbn)< 1e-2
            break
        elseif Thb(N)-Thbn<0
                  Thb2=max(Thb(1),Thb2);
              else
                  Thb1=min(Thb(1),Thb1);
        end
                   Thb(1)=(Thb1+Thb2)/2;


        oo=oo+1;
end


hv0=1000*(2501.689845+1.806916015*(Thb1-273.15)+ ...
                    5.087717*10^-4*(Thb1-273.15)^2-1.1221*10^-5*(Thb1-273.15)^3);
Qp=mhn-mh(1);
Jw=Qp/w/L*3600;
WR=Qp/mhn;
SEC=mhn*Cpl_h*(Thbn-To)/Qp/1000;
GOR=hv0/SEC/1000;
end