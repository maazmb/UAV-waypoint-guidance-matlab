
%variables declarations
%global variables px,py,pz are cartesian positions for calculating horizontal distance using height,lat and lon.
%xea,yea,zea are cartesian positions also considering the ellipticity of earth.

global lata lona xea yea zea pxa pya pza height distancea psia thetaa phia vHa;      %declare in every function
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


lat1 = (33.718151)*pi/180;                 %position1 = islamabad (frm google earth: 33.718151, 73.060547) ,(also frm net 33+45/60,,73+8/60)
lon1 = (73.060547)*pi/180 ;                 %position2 = peshawar   (frm google earth: 34.004299, 71.54483) ,(also frm net 34,,71+35/60)                      
lat2 = (34.004299)*pi/180;                        %GC dist (from internet)   = 145.7671km  AND 146.0606
lon2 = (71.54483)*pi/180;


xe1 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat1)))*cos(lat1)*cos(lon1);
ye1 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat1)))*cos(lat1)*sin(lon1);
ze1 = (r_e*(1-epsi*epsi)/sqrt(1-epsi*epsi*sin(lat1)*sin(lat1)))*sin(lat1);

xe2 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat2)))*cos(lat2)*cos(lon2);
ye2 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat2)))*cos(lat2)*sin(lon2);
ze2= (r_e*(1-epsi*epsi)/sqrt(1-epsi*epsi*sin(lat1)*sin(lat2)))*sin(lat2);

d = sqrt((xe1-xe2)*(xe1-xe2) + (ye1-ye2)*(ye1-ye2) + (ze1-ze2)*(ze1-ze2));
fprintf('total distance =  %f\n',d);     %146.057 km


%%Distance by using altitude,lat and lon

px1 = r_e * cos(lat1) * sin(lon1);
py1 = r_e * sin(lat1);
pz1 = r_e * cos(lat1) * cos(lon1);

px2 = r_e * cos(lat2) * sin(lon2);
py2 = r_e * sin(lat2);
pz2 = r_e * cos(lat2) * cos(lon2);

dh = sqrt((px2-px1)*(px2-px1) + (py2-py1)*(py2-py1) + (pz2-pz1)*(pz2-pz1));
fprintf('total horizontal distance =  %f\n',dh);                                       % 145.930 km

%%{

nbr = 5/deltat;                           %changed from orignal code
lat2t = lat2;                                %(33+49/60)*pi/180;
lon2t = lon2;                                %(72)*pi/180;
ntax = 300;                                  % 1.2*60/deltat                     %ntax = 720

vitax = 0; 
vftax = 100;
psi1 = 282.80837947*pi/180;


acceltraj(lat1,lon1,lat2t,lon2t,ntax,0,vftax,vitax,deltat,r_e,psi1,0,0,0,1.095,0,vitax*cos(psi1),vitax*sin(psi1),0,0,0,0,0,xe1,ye1,ze1); %also passed extra disatance(ntotal) variable
ntotal = ntax;




%......................................................................................................................
%AlTERATION ALGORITHM COMPUTING THE BEARING 
deg=1;
nm=60.02*deg;
m=1852*nm;

%
%NOTE: I have added this part to compute the required heading between two points, so that the heading calculated could be used as psi1.
%NOTE: Rhumb line bearing came out to  be 86 deg...but using it as psi still didnt do any good to results. 
%NOTE: if Great Circle bearing (-86 deg) used as psi, it moved the lat lon in the right direction but very very slowly. 

%great circle computations
distancegc=2*asin(sqrt(func1(lat2t-lat1) + cos(lat1).*cos(lat2t).*func1(lon2t-lon1)));
    distancegc= distancegc*180*m/pi;
    fprintf('great circle distance=  %0.8f\n',distancegc);   % 145.720 km
    
    bearinggc = atan2(sin(lon2t-lon1).*cos(lat2t),cos(lat1).*sin(lat2t) - sin(lat1).*cos(lat2t).*cos(lon2t-lon1))*180/pi;
    fprintf('great circle bearing=  %0.8f\n',bearinggc);
%}
    
    %

%rhumb line
     y=log(tan(lat2t./2+pi/4)./tan(lat1./2+pi/4));
  if (lat2t==lat1)
     v = cos(lat1);
   else 
     v = (lat2t-lat1)/y;
  end   
  if (lon2t > lon1)
      course = mod(180.*(atan2((lon2t-lon1),y))./pi,360);
      totaldistance= 180*m.*(sqrt(v^2*(lon2t-lon1)^2 + (lat2t-lat1)^2))./pi;
   else
      course = mod(180.*(atan2((lon2t-lon1),y))./pi,360);
      totaldistance= 180*m.*(sqrt(v.^2.*(lon1-lon2t).^2 + (lat2t-lat1).^2))./pi;
  end
     fprintf('Rhumb line distance =  %0.8f\n',totaldistance);            %145.721 km
     fprintf('Rhumb line course =  %0.8f\n',course);
     
    %}


%........................................................................................................................


%changed from orignal code
lat2acc1 = lat2;                                %(33+49/60)*pi/180;
lon2acc1 = lon2;                                %(72)*pi/180;
nacc1 = 400;                                 % 1.2*60/deltat                     %ntax = 720
viacc1 = 0; 
vfacc1 = 1250;

%282.80837947  isl-peshawar rhumb line bearing
%acceltraj(lat1,lon1,lat2acc1,lon2acc1,nacc1,0,vfacc1,viacc1,deltat,r_e,psi1,0,0,0,1900,0,viacc1*cos(psi1),viacc1*sin(psi1),0,0,0,0,0,xe1,ye1,ze1);
acceltraj(lata(ntotal),lona(ntotal),lat2t,lon2t,nacc1,ntotal,vfacc1,viacc1,deltat,r_e,trjd.psi,height(ntotal),0,distancea(ntotal),2100,0,trjd.vx,trjd.vy,trjd.theta,trjd.phi,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal)); %also passed extra disatance(ntotal) variable
ntotal = nacc1 + ntotal;
%distancea(ntax)
trjd
%}
%{
figure
subplot(2,2,1)
plot((lona.*180/pi),(lata.*180/pi))
title('longitude vs latitude')
xlabel('longitude')
ylabel('latitude')

subplot(2,2,2)
plot((distancea),(height))
title('distance vs height')
xlabel('distance')
ylabel('height')

subplot(2,2,3)
plot((trjd.lon.*180/pi),(trjd.psi.*180/pi))
title('longitude vs psi')
xlabel('lon')
ylabel('psi')

subplot(2,2,4)
plot((trjd.lon.*180/pi),(trjd.vfin.*180/pi))
title('longitude vs vh')
xlabel('lon')
ylabel('vfin')
%}
%%{
nacc2 = 250;                           %changed from orignal code
lat2acc2 = lat2;                                %(33+49/60)*pi/180;
lon2acc2 = lon2;                                %(72)*pi/180;
                                 % 1.2*60/deltat                     %ntax = 720
viacc2 = 0; 
vfacc2 = 2000;



%acceltraj(lat1,lon1,lat2acc2,lon2acc2,nacc2,0,vfacc2,viacc2,deltat,r_e,psi1,0,0,0,0,0,viacc2*cos(psi1),viacc2*sin(psi1),0,0,0,0,0,xe1,ye1,ze1);
acceltraj(lata(ntotal),lona(ntotal),lat2acc2,lon2acc2,nacc2,ntotal,vfacc2,viacc2,deltat,r_e,trjd.psi,height(ntotal),0,distancea(ntotal),0,0,trjd.vx,trjd.vy,trjd.theta,trjd.phi,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal)); %also passed extra disatance(ntotal) variable
%deacceltraj(lata(ntotal),lona(ntotal),lat2cl,lon2cl,nclimb,ntotal,vfcl,trjd.vfin,deltat,r_e,trjd.psi,height(ntotal),0.3048,distancea(ntotal),200,0,trjd.vx,trjd.vy,trjd.theta,trjd.phi,vfx,vfy,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal))

ntotal = nacc2 + ntotal;
trjd

%{
figure
subplot(2,2,1)
plot((lona.*180/pi),(lata.*180/pi))
title('longitude vs latitude')
xlabel('longitude')
ylabel('latitude')

subplot(2,2,2)
plot((distancea),(height))
title('distance vs height')
xlabel('distance')
ylabel('height')

subplot(2,2,3)
plot((trjd.lon.*180/pi),(trjd.psi.*180/pi))
title('longitude vs psi')
xlabel('lon')
ylabel('psi')

subplot(2,2,4)
plot((trjd.lon.*180/pi),(trjd.vfin.*180/pi))
title('longitude vs vh')
xlabel('lon')
ylabel('vfin')

filename = 'file1.kml';
kmlwrite(filename,(lata*180/pi),(lona*180/pi),height);
%winopen(filename)
%}


%%kml file for google earth

nacc3 = 25;                           %changed from orignal code
lat2acc3 = lat2;                                %(33+49/60)*pi/180;
lon2acc3 = lon2;                                %(72)*pi/180;
                                 % 1.2*60/deltat                     %ntax = 720
viacc3 = 0; 
vfacc3 = -100;



%acceltraj(lat1,lon1,lat2acc2,lon2acc2,nacc2,0,vfacc2,viacc2,deltat,r_e,psi1,0,0,0,0,0,viacc2*cos(psi1),viacc2*sin(psi1),0,0,0,0,0,xe1,ye1,ze1);
acceltraj(lata(ntotal),lona(ntotal),lat2acc3,lon2acc3,nacc3,ntotal,vfacc3,viacc3,deltat,r_e,trjd.psi,height(ntotal),400,distancea(ntotal),0,5000,trjd.vx,trjd.vy,trjd.theta,trjd.phi,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal)); %also passed extra disatance(ntotal) variable


ntotal = nacc3 + ntotal;


nacc4 = 60;                           %changed from orignal code
lat2acc4 = lat2;                                %(33+49/60)*pi/180;
lon2acc4 = lon2;                                %(72)*pi/180;
                                 % 1.2*60/deltat                     %ntax = 720
viacc4 = 0; 
vfacc4 = 100;



%acceltraj(lat1,lon1,lat2acc2,lon2acc2,nacc2,0,vfacc2,viacc2,deltat,r_e,psi1,0,0,0,0,0,viacc2*cos(psi1),viacc2*sin(psi1),0,0,0,0,0,xe1,ye1,ze1);
acceltraj(lata(ntotal),lona(ntotal),lat2acc4,lon2acc4,nacc4,ntotal,vfacc4,viacc4,deltat,r_e,trjd.psi,height(ntotal),0,distancea(ntotal),0,240,trjd.vx,trjd.vy,trjd.theta,trjd.phi,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal)); %also passed extra disatance(ntotal) variable


ntotal = nacc4 + ntotal;



nacc5 = 65;                           %changed from orignal code
lat2acc5 = lat2;                                %(33+49/60)*pi/180;
lon2acc5 = lon2;                                %(72)*pi/180;
                                 % 1.2*60/deltat                     %ntax = 720
viacc5 = 0; 
vfacc5 = -100;



%acceltraj(lat1,lon1,lat2acc2,lon2acc2,nacc2,0,vfacc2,viacc2,deltat,r_e,psi1,0,0,0,0,0,viacc2*cos(psi1),viacc2*sin(psi1),0,0,0,0,0,xe1,ye1,ze1);
acceltraj(lata(ntotal),lona(ntotal),lat2acc5,lon2acc5,nacc5,ntotal,vfacc5,viacc5,deltat,r_e,trjd.psi,height(ntotal),200,distancea(ntotal),0,4000,trjd.vx,trjd.vy,trjd.theta,trjd.phi,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal)); %also passed extra disatance(ntotal) variable


ntotal = nacc5 + ntotal;


nacc6 = 9;                           %changed from orignal code
lat2acc6 = lat2;                                %(33+49/60)*pi/180;
lon2acc6 = lon2;                                %(72)*pi/180;
                                 % 1.2*60/deltat                     %ntax = 720
viacc6 = 0; 
vfacc6 = 0;



%acceltraj(lat1,lon1,lat2acc2,lon2acc2,nacc2,0,vfacc2,viacc2,deltat,r_e,psi1,0,0,0,0,0,viacc2*cos(psi1),viacc2*sin(psi1),0,0,0,0,0,xe1,ye1,ze1);
acceltraj(lata(ntotal),lona(ntotal),lat2acc6,lon2acc6,nacc6,ntotal,vfacc6,viacc6,deltat,r_e,trjd.psi,height(ntotal),390.9,distancea(ntotal),0,300,trjd.vx,trjd.vy,trjd.theta,trjd.phi,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal)); %also passed extra disatance(ntotal) variable


ntotal = nacc6 + ntotal;

nacc7 = 11;                           %changed from orignal code
lat2acc7 = lat2;                                %(33+49/60)*pi/180;
lon2acc7 = lon2;                                %(72)*pi/180;
                                 % 1.2*60/deltat                     %ntax = 720
viacc7 = 0; 
vfacc7 = 0;



%acceltraj(lat1,lon1,lat2acc2,lon2acc2,nacc2,0,vfacc2,viacc2,deltat,r_e,psi1,0,0,0,0,0,viacc2*cos(psi1),viacc2*sin(psi1),0,0,0,0,0,xe1,ye1,ze1);
acceltraj(lata(ntotal),lona(ntotal),lat2acc7,lon2acc7,nacc7,ntotal,vfacc7,viacc7,deltat,r_e,trjd.psi,height(ntotal),0,distancea(ntotal),0,19,trjd.vx,trjd.vy,trjd.theta,trjd.phi,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal)); %also passed extra disatance(ntotal) variable


ntotal = nacc7 + ntotal;







%%inorder to test the acceleration function only ,put breakpoint here
%psi1 = 282.80837947*pi/180;

trjd

figure
plot((lona.*180/pi),(lata.*180/pi))
title('longitude vs latitude')
xlabel('longitude')
ylabel('latitude')

figure
plot((distancea),(height))
title('distance vs height')
xlabel('distance')
ylabel('height')

figure
plot((lona.*180/pi),(psia.*180/pi))
title('longitude vs psia')
xlabel('lon')
ylabel('psi')

figure
plot((lona.*180/pi),(vHa))
title('longitude vs vh')
xlabel('lon')
ylabel('vH')

height(ntotal)


% for simulink
time = 0:1:1119;
gdlata = geoc2geod((lata*180/pi),6378137,'WGS84');
lonasim=[time;(lona*180/pi)]';
latasim=[time;gdlata]';
heightsim=[time;height]';
distanceasim = [time;distancea]';
phiasim=[time;(phia*180/pi)]';
thetaasim=[time;(thetaa*180/pi)]';
psiasim=[time;(psia*180/pi)]';

%phiasims= single(phiasim);
%thetaasims=single(thetaasim);
%psiasims=single(psiasim);
%psiasingle = single(psia*180/pi);


filename = 'file1.kml';
kmlwrite(filename,(lata*180/pi),(lona*180/pi),height);




%{
ndeacc1 = 100;                         
lat2deacc1 = lat2;
lon2deacc1 = lon2;
videacc1=0;                              %QUESTION :: its never used?
vfdeacc1 = 0;
vfx = -50;
vfy = -200;

%deacceltraj(lat1,lon1,lat2,lon2,nclimb,0,vfcl,vicl,deltat,r_e,psi1,0,0,0,0,4000,vicl*cos(psi1),vicl*sin(psi1),0,0,vfx,vfy,0,0,0,xe1,ye1,ze1)
deacceltraj(lata(ntotal),lona(ntotal),lat2deacc1,lon2deacc1,ndeacc1,ntotal,vfdeacc1,trjd.vfin,deltat,r_e,trjd.psi,height(ntotal),0,distancea(ntotal),0,4000,trjd.vx,trjd.vy,trjd.theta,trjd.phi,vfx,vfy,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal))
ntotal = ntotal + ndeacc1;
%}


nturn1=10
psi2 = 110*pi/180;
height(ntotal)=0;
turntraj(lata(ntotal),lona(ntotal),trjd.psi,psi2,nturn1,ntotal,trjd.vfin,trjd.vfin,deltat,r_e,height(ntotal),trjd.fb,trjd.vx,trjd.vy,trjd.theta,trjd.phi,1000,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal))
ntotal = ntotal+nturn1+nbr;

%{
%%ANALYSIS::  
1) Even though the computed lat and lons are continous but they are not in  the right  direction. If latitude has to be decreased 
 in acceleration..the output gives increasing latitudes,,etc

2) I feel this is maybe because the acceleration function isnt using the    final lat and lon passed (lat2t,lon2t) anywhere
 on its updates of variables.

3) I thought since the updates depend on heading (psi) so i should try  to     find a specific psi to make it move in the right direction.
 so rather than putting random value of psi i tested using bearing from
 rhumb line and great circle formulas.
 


%}

