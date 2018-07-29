function [  ] = turntraj(lat1,lon1,psi1,psi2,n1,n2,vf,vi,deltat,r_e,hei,fb,vx,vy,th1,ph1,turnradius,vxe,vye,vze,xei,yei,zei)
global trjd;
global lata lona xea yea zea height;
global epsi omega;

negtr = 1;
if (psi2 < psi1)
    negtr = -1;
end

a = (vf-vi)/(n1);         % QUESTION :: not(n1*deltat)
%fb = a;
wx = 0;
wy = 0;
wz = 0;
g = 9.8;
R = 10000;
fb = vi*vi/R;

psi(1) = psi1;
theta(1) = th1;
phi(1) = ph1;

vN(1) = vx;
vE(1) = vy;
%vD(1) = 0;

vxep = vxe;
vyep = vye;
vzep = vze;

lat(1) = lat1;
lon(1) = lon1;
xe(1) = xei;
ye(1) = yei;
ze(1) = zei;


%height(n2+1) = h1;
lata(n2+1)= lat1;
display(lata(n2+1))
lona(n2+1)= lon1;
display(lona(n2+1))
xea(n2+1) = xei;
yea(n2+1) = yei;
zea(n2+1) = zei;



for k=2:(n1) 
    
      wz = negtr*3*pi/180;  %3 deg per sec
    
      phi(k)=phi(k-1)+((wy*sin(phi(k-1))+wz*cos(phi(k-1)))*tan(theta(k-1))+wx)*deltat;
      theta(k)=theta(k-1)+(wy*cos(phi(k-1))-wz*sin(phi(k-1)))*deltat;
      psi(k)=psi(k-1)+((wy*sin(phi(k-1))+wz*cos(phi(k-1)))/(cos(theta(k-1))))*deltat;

      fbx = 0;
      fby = wz*wz*turnradius;
      
      fN= -fby*sin(psi(k));
      fE=  fby*cos(psi(k));
 
      %vN(k)=vN(k-1)+(fN-2*Omega*vE(k-1)*sin(lata(k))+(vN(k-1)*vD(k-1)-vE(k-1)*vE(k-1)*tan(lata(k-1)))/(r_e))*deltat;
      %vE(k)=vE(k-1)+(fE-2*Omega*(vN(k-1)*sin(lata(k))+vD(k-1)*cos(lata(k)))+vE(k-1)*(vD(k-1)+vN(k-1)*tan(lata(k-1)))/(r_e))*deltat;
      %vN(k)=vN(k-1)+(fN-2*Omega*vE(k-1)*sin(lat(k))+(-vE(k-1)*vE(k-1)*tan(lat(k-1)))/(r_e))*deltat;
      %vE(k)=vE(k-1)+(fE-2*Omega*(vN(k-1)*sin(lat(k)))+vE(k-1)*(vN(k-1)*tan(lat(k-1)))/(r_e))*deltat;
    
     
     vN(k) = vN(k-1) + fbx*deltat;
     vE(k) = vE(k-1) + fby*deltat;
     
     %lat(k)= lat(k-1)+vN(k)/*/(r_e)*/*deltat;
     %lon(k)= lon(k-1)+vE(k)/*/(r_e*cos(lat(k)))*/*deltat;
     %vN
     %l=vN /r_e*deltat
     lat(k)= lat(k-1)+vN(k) /r_e*deltat;
     lon(k)= lon(k-1)+vE(k)/(r_e*cos(lat(k)))*deltat;
     %h(k)=h(k-1)-vD(k)*deltat;

     
     vxe=vxep+deltat*(-fN*sin(lat(k))*sin(lon(k))-fE*sin(lat(k))*sin(lon(k))+2*omega*vyep+g*cos(lat(k)));
     vye=vyep+deltat*(-fN*sin(lon(k))+fE*cos(lon(k))-2*omega*vxep);
     vze=vzep+deltat*(-fN*cos(lat(k))*cos(lon(k))-fE*cos(lat(k))*cos(lon(k))-g*sin(lat(k)) );
     vxep=vxe;
     vyep=vye;
     vzep=vze;
     
     vH=sqrt(vN(k)*vN(k)+vE(k)*vE(k));
     
     %vD(k) = vD(k-1);     
    
     %xe(k-1)
    %vxe
     xe(k)=xe(k-1)+vxe*deltat;
     ye(k)=ye(k-1)+vye*deltat;
     ze(k)=ze(k-1)+vze*deltat; 
     xe(k)=(r_e/sqrt(1-epsi*epsi*sin(lat(k))*sin(lat(k))) )*cos(lat(k))*cos(lon(k));
     ye(k)=(r_e/sqrt(1-epsi*epsi*sin(lat(k))*sin(lat(k))))*cos(lat(k))*sin(lon(k));
     ze(k)=(r_e*(1-epsi*epsi)/sqrt(1-epsi*epsi*sin(lat(k))*sin(lat(k))))*sin(lat(k));
     

     lata(n2+k)=lat(k);
     lona(n2+k)=lon(k);
     %display(k);
     %display(lata(n2+k));
     %display(lona(n2+k));
     
     
     
     xea(n2+k)=xe(k);
     yea(n2+k)=ye(k);
     zea(n2+k)=ze(k);
     %height(n2+k)=hei;
    
    
end


%2nd for loop......................

fby = (0-vE(n1))/5;
fbx = 0;
nbr = 5/deltat;

for k=(n1+1):(n1+nbr)
    
      phi(k)=phi(k-1);
      theta(k)=theta(k-1);
      psi(k)=psi(k-1);
      
      %fby = wz*wz*turnradius;
      
      fN= -fby*sin(psi(k));
      fE=  fby*cos(psi(k));

      vN(k) = vN(k-1) + fbx*deltat;
      vE(k) = vE(k-1) + fby*deltat;

     lat(k)= lat(k-1)+vN(k) /r_e*deltat;
     lon(k)= lon(k-1)+vE(k)/(r_e*cos(lat(k)))*deltat;
     
     vxe=vxep+deltat*(-fN*sin(lat(k))*sin(lon(k))-fE*sin(lat(k))*sin(lon(k))+2*omega*vyep+g*cos(lat(k)));
     vye=vyep+deltat*(-fN*sin(lon(k))+fE*cos(lon(k))-2*omega*vxep);
     vze=vzep+deltat*(-fN*cos(lat(k))*cos(lon(k))-fE*cos(lat(k))*cos(lon(k))-g*sin(lat(k)) );
     vxep=vxe;
     vyep=vye;
     vzep=vze;
     
      vH=sqrt(vN(k)*vN(k)+vE(k)*vE(k));
     
     %vD(k) = vD(k-1);     
    
     %xe(k-1)
    %vxe
     xe(k)=xe(k-1)+vxe*deltat;
     ye(k)=ye(k-1)+vye*deltat;
     ze(k)=ze(k-1)+vze*deltat; 
     xe(k)=(r_e/sqrt(1-epsi*epsi*sin(lat(k))*sin(lat(k))) )*cos(lat(k))*cos(lon(k));
     ye(k)=(r_e/sqrt(1-epsi*epsi*sin(lat(k))*sin(lat(k))))*cos(lat(k))*sin(lon(k));
     ze(k)=(r_e*(1-epsi*epsi)/sqrt(1-epsi*epsi*sin(lat(k))*sin(lat(k))))*sin(lat(k));
     

     lata(n2+k)=lat(k);
     lona(n2+k)=lon(k);
     %display(k);
     %display(lata(n2+k));
     %display(lona(n2+k));
     
     
     
     xea(n2+k)=xe(k);
     yea(n2+k)=ye(k);
     zea(n2+k)=ze(k);
     %height(n2+k)=hei;
    
end



%lat = lat1+1;
%lon = lon1+1;
%vH=9;
%psi=psi0+1;
%fb=0;
%vN=0;
%vE=0;
%phi=0;  
%theta=0;
%vxe=vup1+1;
%vye=0;
%vze=0;

display(lat);
display(lon);
trjd.lat = lat;
trjd.lon = lon;
trjd.vfin = vH;
trjd.psi = psi(n1+nbr);
trjd.fb = fbx;
trjd.vx = vN(n1+nbr);
trjd.vy = vE(n1+nbr);
trjd.phi = phi(n1+nbr);
trjd.theta = theta(n1+nbr);
trjd.vhor = vH;
trjd.vex = vxe;
trjd.vey = vye;
trjd.vez = vze;





end

