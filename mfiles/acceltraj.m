function [] = acceltraj(lat1,lon1,lat2,lon2,n1,n2,vf,vi,deltat,r_e,psi1,h1,h2,vup1,vup2,vx,vy,th1,ph1,vxe,vye,vze,xei,yei,zei)

global trjd;
global lata lona xea yea zea height;
global epsi omega;



pox1 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat1)))*cos(lat1)*cos(lon1);
poy1 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat1)))*cos(lat1)*sin(lon1);
poz1= (r_e*(1-epsi*epsi)/sqrt(1-epsi*epsi*sin(lat1)*sin(lat1)))*sin(lat1);

pox2 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat2)))*cos(lat2)*cos(lon2);
poy2 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat2)))*cos(lat2)*sin(lon2);
poz2= (r_e*(1-epsi*epsi)/sqrt(1-epsi*epsi*sin(lat1)*sin(lat2)))*sin(lat2);

d = sqrt((pox1-pox2)*(pox1-pox2) + (poy1-poy2)*(poy1-poy2));
fprintf('distance =  %f\n',d);

a = (vf-vi)/(n1*deltat)
fb = a;
fbx = fb
fby = 0;
fD = -9.8 + (vup2-vup1)/n1;
g = 9.8;

wx = 0;
wy = 0;
wz = 0;

vxep = vxe;
vyep = vye;
vzep = vze;

psi(1) = psi1;
theta(1) = th1;
phi(1) = ph1;

%vNp = vx;
%vEp = vy;
vN(1) = vx;
vE(1) = vy;
vD(1) = vup1; 
vD(1) = (h2-h1/n1*deltat);
h(1) = h1;

lat(1) = lat1;
lon(1) = lon1;
xe(1) = xei;
ye(1) = yei;
ze(1) = zei;
height(n2+1) = h1;
lata(n2+1)= lat1;
display(lata(n2+1))
lona(n2+1)= lon1;
display(lona(n2+1))
xea(n2+1) = xei;
yea(n2+1) = yei;
zea(n2+1) = zei;


for k=2:(n1) 
    
      phi(k)=phi(k-1)+((wy*sin(phi(k-1))+wz*cos(phi(k-1)))*tan(theta(k-1))+wx)*deltat;
      theta(k)=theta(k-1)+(wy*cos(phi(k-1))-wz*sin(phi(k-1)))*deltat;
      psi(k)=psi(k-1)+((wy*sin(phi(k-1))+wz*cos(phi(k-1)))/(cos(theta(k-1))))*deltat;

     
      fN=fb*cos(theta(k))*cos(psi(k));
      fE=fb*cos(theta(k))*sin(psi(k));
 
      %% QUESTION :: the orignal formula doesnt contain (1) (r_e+h) and it doesnt considers (2)the effect of gravity. (3) and sin(lata(k)) while tan(lata(k-1))? why not same angle (k-1)?? 
      %% NOTE : i have assumed (k-1)th of all lat and lon. And also inserted height which seems to have no effect on results anyways.
      
      %vN(k)=vN(k-1)+(fN-2*Omega*vE(k-1)*sin(lata(k))+(vN(k-1)*vD(k-1)-vE(k-1)*vE(k-1)*tan(lata(k-1)))/(r_e))*deltat;
      %vE(k)=vE(k-1)+(fE-2*Omega*(vN(k-1)*sin(lata(k))+vD(k-1)*cos(lata(k)))+vE(k-1)*(vD(k-1)+vN(k-1)*tan(lata(k-1)))/(r_e))*deltat;
      %vN(k)=vN(k-1)+(fN-2*Omega*vE(k-1)*sin(lat(k))+(vN(k-1)*vD(k-1)-vE(k-1)*vE(k-1)*tan(lat(k-1)))/(r_e))*deltat;
      %vE(k)=vE(k-1)+(fE-2*Omega*(vN(k-1)*sin(lat(k))+vD(k-1)*cos(lat(k)))+vE(k-1)*(vD(k-1)+vN(k-1)*tan(lat(k-1)))/(r_e))*deltat;
      %vD(k)=vD(k-1)+(fD-2*Omega*vE(k)*cos(lat(k))-(vE(k)*vE(k)+vN(k)*vN(k))/(r_e+h(k))+9.8)*deltat;  
      
      vN(k)=vN(k-1)+(fN-2*omega*vE(k-1)*sin(lat(k-1))+(vN(k-1)*vD(k-1)-vE(k-1)*vE(k-1)*tan(lat(k-1)))/(r_e+h(k-1) ))*deltat;
      vE(k)=vE(k-1)+(fE-2*omega*(vN(k-1)*sin(lat(k-1))+vD(k-1)*cos(lat(k-1)))+vE(k-1)*(vD(k-1)+vN(k-1)*tan(lat(k-1)))/(r_e+h(k-1)))*deltat;
      vD(k)=vD(k-1)+(fD-2*omega*vE(k-1)*cos(lat(k-1))-(vE(k-1)*vE(k-1)+vN(k-1)*vN(k-1))/(r_e+h(k-1))+9.8)*deltat;
    
      
     %% NOTE: i have not used the following four formulas and used instead the above ones.Both give the same result. 
     
     %vN=vNp+fbx*cos(psi(k))*deltat;
     %vE=vEp+fbx*sin(psi(k))*deltat;
     %vEp=vE;
     %vNp=vN;
     vH=sqrt(vN(k)*vN(k)+vE(k)*vE(k));
     
     
     %lat(k)= lat(k-1)+vN(k)/*/(r_e)*/*deltat;
     %lon(k)= lon(k-1)+vE(k)/*/(r_e*cos(lat(k)))*/*deltat;
     %vN
     %l=vN /r_e*deltat
     lat(k)= lat(k-1)+vN(k) /r_e*deltat;
     lon(k)= lon(k-1)+vE(k)/(r_e*cos(lat(k)))*deltat;
     h(k)=h(k-1)-vD(k)*deltat;

     
     vxe=vxep+deltat*(-fN*sin(lat(k))*sin(lon(k))-fE*sin(lat(k))*sin(lon(k))+2*omega*vyep+g*cos(lat(k)));
     vye=vyep+deltat*(-fN*sin(lon(k))+fE*cos(lon(k))-2*omega*vxep);
     vze=vzep+deltat*(-fN*cos(lat(k))*cos(lon(k))-fE*cos(lat(k))*cos(lon(k))-g*sin(lat(k)) );
     vxep=vxe;
     vyep=vye;
     vzep=vze;
    %xe(k-1)
    %vxe
     xe(k)=xe(k-1)+vxe*deltat;
     ye(k)=ye(k-1)+vye*deltat;
     ze(k)=ze(k-1)+vze*deltat; 

     

     lata(n2+k)=lat(k);
     lona(n2+k)=lon(k);
     %display(k);
     %display(lata(n2+k));
     %display(lona(n2+k));
     
     xe(k)=(r_e/sqrt(1-epsi*epsi*sin(lat(k))*sin(lat(k))) )*cos(lat(k))*cos(lon(k));
     ye(k)=(r_e/sqrt(1-epsi*epsi*sin(lat(k))*sin(lat(k))))*cos(lat(k))*sin(lon(k));
     ze(k)=(r_e*(1-epsi*epsi)/sqrt(1-epsi*epsi*sin(lat(k))*sin(lat(k))))*sin(lat(k)); 
     
     xea(n2+k)=xe(k);
     yea(n2+k)=ye(k);
     zea(n2+k)=ze(k);
     height(n2+k)=h(k);
    
    
end

%{
lat = lat1+1;
lon = lon1+1;
vH=9;
psi=psi0+1;
fb=0;
vN=0;
vE=0;
phi=0;  
theta=0;
vxe=vup1+1;
vye=0;
vze=0;
%}

fprintf('initial latitude passed =  %0.8f\n',lat1);
fprintf('initial longitude passed =  %0.8f\n',lon1);
fprintf('final latitude passed =  %0.8f\n',lat2);
fprintf('final longitude passed =  %0.8f\n',lon2);
fprintf('final latitude =  %0.8f\n',lata(n2+k));
fprintf('final longitude =  %0.8f\n',lona(n2+k));
fprintf('final height =  %0.8f\n',height(n2+k));

trjd.lat = lat;
trjd.lon = lon;
trjd.vfin = vH;
trjd.psi = psi(n1);
trjd.fb = fb;
trjd.vx = vN(n1);
trjd.vy = vE(n1);
trjd.phi = phi(n1);
trjd.theta = theta(n1);
trjd.vhor = vH;
trjd.vex = vxe;
trjd.vey = vye;
trjd.vez = vze;


end

