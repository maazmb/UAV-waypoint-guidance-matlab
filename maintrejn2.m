%flight trejactory from KSFO Airport to KSQL San Carlos Airport 
%trejactory 1

%variables declarations
%global variables px,py,pz are cartesian positions for calculating horizontal distance using height,lat and lon.
%xea,yea,zea are cartesian positions also considering the ellipticity of earth.

global lata lona xea yea zea pxa pya pza height distancea psia thetaa phia vHa vxea vyea vzea  vNa vEa vDa vNED;      %declare in every function
global trjd;
global epsi; 
global omega;                             % = 7.2922116*10^(-5)



f = 298.257223563;                        %ellipticity
omega = 7.2922116*10^(-5);
r_e = 6378137;
r_p = r_e*(1-1/f);
epsi = sqrt(1-(r_p*r_p)/(r_e*r_e));
%epsi = 0;
deltat = 0.1;
g = 9.8;
FL = 280;


lat1 = (37.6287)*pi/180;                   %position1 = KSFO Sans Francisco Airport runway 10l begining 
lon1 = (-122.393)*pi/180 ;                 %position2 = KSFO Sans Francisco Airport runway 10l final                      
lat2 = (37.616452)*pi/180;                        %GC dist = 2.8967 km  AND rhumb line bearing= 118.0334 deg
lon2 = (-122.363958)*pi/180;


xe1 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat1)))*cos(lat1)*cos(lon1);
ye1 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat1)))*cos(lat1)*sin(lon1);
ze1 = (r_e*(1-epsi*epsi)/sqrt(1-epsi*epsi*sin(lat1)*sin(lat1)))*sin(lat1);

xe2 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat2)))*cos(lat2)*cos(lon2);
ye2 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat2)))*cos(lat2)*sin(lon2);
ze2= (r_e*(1-epsi*epsi)/sqrt(1-epsi*epsi*sin(lat1)*sin(lat2)))*sin(lat2);

d = sqrt((xe1-xe2)*(xe1-xe2) + (ye1-ye2)*(ye1-ye2) + (ze1-ze2)*(ze1-ze2));
fprintf('total distance =  %f\n',d);     % km


%%Distance by using altitude,lat and lon

px1 = r_e * cos(lat1) * sin(lon1);
py1 = r_e * sin(lat1);
pz1 = r_e * cos(lat1) * cos(lon1);

px2 = r_e * cos(lat2) * sin(lon2);
py2 = r_e * sin(lat2);
pz2 = r_e * cos(lat2) * cos(lon2);

dh = sqrt((px2-px1)*(px2-px1) + (py2-py1)*(py2-py1) + (pz2-pz1)*(pz2-pz1));
fprintf('total horizontal distance =  %f\n',dh);                                       %  km

%%{

nbr = 5/deltat;                           
lat2t = lat2;                                
lon2t = lon2;                                
ntax = 20;                                  % 1.2*60/deltat                     %ntax = 720
%1400 30
vitax = 0; 
vftax = 30;
psi1 = 118.0334*pi/180;    %rhumb line heading
%118.0334 2.286


acceltraj(lat1,lon1,lat2t,lon2t,ntax,0,vftax,vitax,deltat,r_e,psi1,2.345696,0,0,-0.12,0,vitax*cos(psi1),vitax*sin(psi1),0,0,0,0,0,xe1,ye1,ze1); %also passed extra disatance(ntotal) variable
ntotal = ntax;


ntax1 = 670;                                  % 1.2*60/deltat                     %ntax = 720
%1400 30
vitax1 = 30; 
vftax1 = 30;

acceltraj(lata(ntotal),lona(ntotal),lat2t,lon2t,ntax1,ntotal,vftax1,vitax1,deltat,r_e,trjd.psi,height(ntotal),0,distancea(ntotal),-1.86,0,trjd.vx,trjd.vy,trjd.theta,trjd.phi,vxea(ntotal),vyea(ntotal),vzea(ntotal),xea(ntotal),yea(ntotal),zea(ntotal)); %also passed extra disatance(ntotal) variable
ntotal = ntax1 + ntotal ;      
   
figure
plot((distancea),(height))
title('distance vs height')
xlabel('distance (m)')
ylabel('height (m)')




ninitclimb = 220;                                 % 1.2*60/deltat                     %ntax = 720
viinitclimb = 30; 
vfinitclimb = 30;

%acceltraj(lat1,lon1,lat2acc1,lon2acc1,nacc1,0,vfacc1,viacc1,deltat,r_e,psi1,0,0,0,1900,0,viacc1*cos(psi1),viacc1*sin(psi1),0,0,0,0,0,xe1,ye1,ze1);
acceltraj(lata(ntotal),lona(ntotal),lat2t,lon2t,ninitclimb,ntotal,vfinitclimb,viinitclimb,deltat,r_e,trjd.psi,height(ntotal),0,distancea(ntotal),880,0,trjd.vx,trjd.vy,trjd.theta,trjd.phi,vxea(ntotal),vyea(ntotal),vzea(ntotal),xea(ntotal),yea(ntotal),zea(ntotal)); %also passed extra disatance(ntotal) variable
ntotal = ninitclimb + ntotal;







nclimb = 550;                         
viclimb = 30; 
vfclimb = 30;



%acceltraj(lat1,lon1,lat2acc2,lon2acc2,nacc2,0,vfacc2,viacc2,deltat,r_e,psi1,0,0,0,0,0,viacc2*cos(psi1),viacc2*sin(psi1),0,0,0,0,0,xe1,ye1,ze1);
acceltraj(lata(ntotal),lona(ntotal),lat2t,lon2t,nclimb,ntotal,vfclimb,viclimb,deltat,r_e,trjd.psi,height(ntotal),0,distancea(ntotal),400,0,trjd.vx,trjd.vy,trjd.theta,trjd.phi,vxea(ntotal),vyea(ntotal),vzea(ntotal),xea(ntotal),yea(ntotal),zea(ntotal)); %also passed extra disatance(ntotal) variable
%deacceltraj(lata(ntotal),lona(ntotal),lat2cl,lon2cl,nclimb,ntotal,vfcl,trjd.vfin,deltat,r_e,trjd.psi,height(ntotal),0.3048,distancea(ntotal),200,0,trjd.vx,trjd.vy,trjd.theta,trjd.phi,vfx,vfy,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal))

ntotal = nclimb + ntotal;
display('after climb')
display(height(ntotal))
display(lata(ntotal)*180/pi)
display(lona(ntotal)*180/pi)
display(distancea(ntotal))





%{
nturn = 250;      
psiturn = 140.739874*pi/180;
viturn = 30; 
vfturn = 30;
turnradius = 1000;

turntraj(lata(ntotal),lona(ntotal),psi1,psiturn,nturn,ntotal,vfturn,viturn,deltat,r_e,height(ntotal),distancea(ntotal),trjd.fb,trjd.vx,trjd.vy,trjd.theta,trjd.phi,turnradius,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal));

ntotal = nturn + ntotal;
%}


ncruise = 3000;      
psi2 = 140.739874*pi/180;
%IF TURN SO 138.454430......TOTAL HEADING WITHOUT TURN 140.739874
vicruise = 30; 
vfcruise = 30;

acceltraj(lata(ntotal),lona(ntotal),lat2t,lon2t,ncruise,ntotal,vfcruise,vicruise,deltat,r_e,psi2,height(ntotal),0,distancea(ntotal),0,15,vicruise*cos(psi2),vicruise*sin(psi2),trjd.theta,trjd.phi,vxea(ntotal),vyea(ntotal),vzea(ntotal),xea(ntotal),yea(ntotal),zea(ntotal)); %also passed extra disatance(ntotal) variable

ntotal = ncruise + ntotal;

display('after cruise')
display(height(ntotal))
display(lata(ntotal)*180/pi)
display(lona(ntotal)*180/pi)
display(distancea(ntotal))


ndesnt = 190;
videsnt =30;
vfdesnt = 30;

acceltraj(lata(ntotal),lona(ntotal),lat2t,lon2t,ndesnt,ntotal,vfdesnt,videsnt,deltat,r_e,trjd.psi,height(ntotal),22,distancea(ntotal),0,0,trjd.vx,trjd.vy,trjd.theta,trjd.phi,vxea(ntotal),vyea(ntotal),vzea(ntotal),xea(ntotal),yea(ntotal),zea(ntotal));
ntotal = ndesnt + ntotal;

display('after desnt')
display(height(ntotal))
display(lata(ntotal)*180/pi)
display(lona(ntotal)*180/pi)
display(distancea(ntotal))


niapp =660;
viiapp =30;
vfiapp = 30;

acceltraj(lata(ntotal),lona(ntotal),lat2t,lon2t,niapp,ntotal,vfiapp,viiapp,deltat,r_e,trjd.psi,height(ntotal),0.5,distancea(ntotal),0,1,trjd.vx,trjd.vy,trjd.theta,trjd.phi,vxea(ntotal),vyea(ntotal),vzea(ntotal),xea(ntotal),yea(ntotal),zea(ntotal));
ntotal = niapp + ntotal;


display('after iinitapp')

display(height(ntotal))
display(lata(ntotal)*180/pi)
display(lona(ntotal)*180/pi)
display(distancea(ntotal))
%for simulink test part delete afterward
%{

time = 0:1:29;
%gdlata = geoc2geod((lata*180/pi),6378137,'WGS84');
lonasim=[time;(lona)]';
latasim=[time;(lata)]';
heightsim=[time;height]';
distanceasim = [time;distancea]';
phiasim=[time;(phia)]';
thetaasim=[time;(thetaa)]';
psiasim=[time;(psia)]';

%phiasims= single(phiasim);
%thetaasims=single(thetaasim);
%psiasims=single(psiasim);
%psiasingle = single(psia*180/pi)
%}






nfapp =680;
vifapp =30;
vffapp = 30;
psi3= 141.581398*pi/180;


acceltraj(lata(ntotal),lona(ntotal),lat2t,lon2t,nfapp,ntotal,vffapp,vifapp,deltat,r_e,psi3,height(ntotal),23.3,distancea(ntotal),0,1,vifapp*cos(psi3),vifapp*sin(psi3),trjd.theta,trjd.phi,vxea(ntotal),vyea(ntotal),vzea(ntotal),xea(ntotal),yea(ntotal),zea(ntotal));
ntotal = nfapp + ntotal;


display('after finaltapp')

display(height(ntotal))
display(lata(ntotal)*180/pi)
display(lona(ntotal)*180/pi)
display(distancea(ntotal))







%filename = 'file1.kml';
%kmlwrite(filename,(lata*180/pi),(lona*180/pi),height);



nland =40;
viland =30;
vfland = 25;
psi4= 159.748655*pi/180;


acceltraj(lata(ntotal),lona(ntotal),lat2t,lon2t,nland,ntotal,vfland,viland,deltat,r_e,psi4,height(ntotal),1,distancea(ntotal),0,0,viland*cos(psi4),viland*sin(psi4),trjd.theta,trjd.phi,vxea(ntotal),vyea(ntotal),vzea(ntotal),xea(ntotal),yea(ntotal),zea(ntotal));
ntotal = nland + ntotal;


%{
nland = 120;
viland= 30;
vfland = 30;
vfx = 25;
vfy = 20;
psi4 = 159.748655*pi/180;
deacceltraj(lata(ntotal),lona(ntotal),lat2t,lon2t,nland,ntotal,vfland,viland,deltat,r_e,psi4,height(ntotal),0.1,distancea(ntotal),0,0,viland*cos(psi4),viland*sin(psi4),trjd.theta,trjd.phi,vfx,vfy,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal));

ntotal = ntotal + nland;

%}
display('after landing')

display(height(ntotal))
display(lata(ntotal)*180/pi)
display(lona(ntotal)*180/pi)
display(distancea(ntotal))
display(vHa(ntotal))







ntax2 = 420;
vitax2= 25.1311;
vftax2 = 0;
vfx2 = 0;
vfy2 = 0;
psi5 = 138.89777*pi/180;
deacceltraj(lata(ntotal),lona(ntotal),lat2t,lon2t,ntax2,ntotal,vftax2,vitax2,deltat,r_e,psi5,height(ntotal),0,distancea(ntotal),-0.75,0,vitax2*cos(psi5),vitax2*sin(psi5),trjd.theta,trjd.phi,vfx2,vfy2,vxea(ntotal),vyea(ntotal),vzea(ntotal),xea(ntotal),yea(ntotal),zea(ntotal));

ntotal = ntotal + ntax2;

  




display('after  taxi 2')

display(height(ntotal))
display(lata(ntotal)*180/pi)
display(lona(ntotal)*180/pi)
display(distancea(ntotal))
display(vHa(ntotal))



figure
plot((lona.*180/pi),(lata.*180/pi))
title('longitude vs latitude')
xlabel('longitude (deg)')
ylabel('latitude (deg)')

figure
plot((distancea),(height))
title('distance vs height')
xlabel('distance (m)')
ylabel('height (m)')

figure
plot((lona.*180/pi),(psia.*180/pi))
title('longitude vs psia')
xlabel('lon (deg)')
ylabel('psi (deg)')



figure
plot((lona.*180/pi),(vHa))
title('longitude vs vh')
xlabel('lon (deg)')
ylabel('vH (m/s)')


figure
plot((distancea),(psia.*180/pi))
title('distance vs heading')
xlabel('distance (m)')
ylabel('psi (deg)')





%}





%}
%....................................................................................................................


%%inorder to test the acceleration function only ,put breakpoint here
%psi1 = 282.80837947*pi/180;

lona = lona;
lata = lata;

vhb = sqrt(vNa.*vNa + vEa.*vEa + vDa.*vDa);
vNED = [vNa;vEa;vDa];

i=1;
dt= 1:i:6450;
phidot=diff(phia);
thetadot=diff(thetaa);
psidot=diff(psia);
phidot(6450)=0;
thetadot(6450)=0;
psidot(6450)=0;


%for velocity in body frame
for n=1:6450 

%t{3,3,6450} = [];    
cn2b(:,:,n) = [cos(psia(n)).*cos(thetaa(n))   sin(psia(n)).*cos(thetaa(n))   -sin(thetaa(n)); sin(thetaa(n)).*cos(psia(n)).*sin(phia(n))-cos(phia(n)).*sin(psia(n))   sin(phia(n)).*sin(psia(n)).*sin(thetaa(n))+cos(psia(n)).*cos(phia(n))   sin(phia(n)).*cos(thetaa(n));  cos(phia(n)).*cos(psia(n)).*sin(thetaa(n))+sin(phia(n)).*sin(psia(n))   sin(psia(n)).*sin(thetaa(n)).*cos(phia(n))-sin(phia(n)).*cos(psia(n))   cos(phia(n)).*cos(thetaa(n))];
%t(:,:,n) = [cos(psia(n))   sin(psia(n))   -sin(thetaa(n)); sin(thetaa(n))   sin(phia(n))   sin(phia(n));  cos(phia(n))  sin(psia(n))   cos(phia(n))];
cb2n(:,:,n) = transpose(cn2b(:,:,n));
vB(:,n) = cn2b(:,:,n)*vNED(:,n);
c1(:,:,n) = [cos(psia(n)) sin(psia(n)) 0; -sin(psia(n)) cos(psia(n)) 0; 0 0 1];
c2(:,:,n) = [cos(thetaa(n)) 0 -sin(thetaa(n));0 1 0; sin(thetaa(n)) 0  cos(thetaa(n))];
c3(:,:,n) = [1 0 0 ;0 cos(phia(n)) sin(phia(n)) ; 0 -sin(phia(n)) cos(phia(n))];
pqr(:,n)= [phidot(n);0;0] + c3(:,:,n)*[0;thetadot(n);0] + c3(:,:,n)*c2(:,:,n)*[0;0;psidot(n)];

end

%i = sqrt((l(2,:)).^2 + (l(1,:)).^2)


vA = sqrt((vB(1,:)).^2 + (vB(2,:)).^2 + (vB(3,:)).^2); 
vAa = sqrt((vB(1,:)).^2 + (vB(2,:)).^2); 




%{
plot((lona.*180/pi),(vA));
figure
plot((lona.*180/pi),(vAa));
figure
plot((lona.*180/pi),(vhb));
figure
plot((lona.*180/pi),(vHa));
%}

%angular rates, psidot on z=r=heading, p = roll=x axis=phi

%%
% RUN THIS PART ONLY AFTER RUNNING SIMULINK MODEL OF DRYDEN IMPLEMENTATION
%comment it while running the code for first time...uncomment it after
%running simulation to run it for second time. 
%{
vBout = transpose(vBout);
pqrout = transpose(pqrout);

for j = 1:6450
    %vB(:,n) = cn2b(:,:,n)*vNED(:,n);
vNEDout(:,j) = cb2n(:,:,j)*vBout(:,j);
end

vNout(:,:) = vNEDout(1,:);
vEout(:,:) = vNEDout(2,:);
vDout(:,:) = vNEDout(3,:);
pout(:,:) = pqrout(1,:);
qout(:,:) = pqrout(2,:);
rout(:,:) = pqrout(3,:);


for z=2:6450 

lata(z)= lata(z-1)+vNout(z) /r_e*deltat;
lona(z)= lona(z-1)+vEout(z)/(r_e*cos(lata(z)))*deltat;
height(z)=height(z-1)-vDout(z)*deltat;
pxa(z) = (r_e+height(z)) * cos(lata(z)) * sin(lona(z));
pya(z) = (r_e+height(z)) * sin(lata(z));
pza(z) = (r_e+height(z)) * cos(lata(z)) * cos(lona(z));
distancea(z) = sqrt((pxa(z)-pxa(1))*(pxa(z)-pxa(1)) + (pya(z)-pya(1))*(pya(z)-pya(1)) + (pza(z)-pza(1))*(pza(z)-pza(1)));

phia(z)=phia(z-1)+((qout(z)*sin(phia(z-1))+rout(z)*cos(phia(z-1)))*tan(thetaa(z-1))+pout(z))*deltat;
thetaa(z)=thetaa(z-1)+(qout(z)*cos(phia(z-1))-rout(z)*sin(phia(z-1)))*deltat;
psia(z)=psi(z-1)+((qout(z)*sin(phia(z-1))+rout(z)*cos(phia(z-1)))/(cos(thetaa(z-1))))*deltat;



end
%angle limiting
psia = wrapTo2Pi(psia);
psia = wrapTo2Pi(psia);
psia = wrapTo2Pi(psia);

%}
%%

% for simulink
time = 0:1:6449;
%gdlata = geoc2geod((lata),6378137,'WGS84');
lonasim=[time;(lona)]';
latasim=[time;lata]';
heightsim=[time;height]';
distanceasim = [time;distancea]';
phiasim=[time;(phia)]';
thetaasim=[time;(thetaa)]';
psiasim=[time;(psia)]';
vHasim = [time;(vHa)]';
rotmatrixsim = [time;psia;thetaa;phia]';
vAsim = [time;vA]';
vBsim = [time;vB]'
pqrsim = [time;pqr]';

%phiasims= single(phiasim);
%thetaasims=single(thetaasim);
%psiasims=single(psiasim);
%psiasingle = single(psia*180/pi)


%filename = 'file1.kml';
%kmlwrite(filename,(lata*180/pi),(lona*180/pi),height);

%Aircraft model: cub














