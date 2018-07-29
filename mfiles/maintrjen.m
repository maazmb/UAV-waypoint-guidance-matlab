
%variables declarations

global lata lona xea yea zea height;      %declare in every function
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


lat1 = (33+45/60)*pi/180;                 %position1 = islamabad
lon1 = (73+8/60)*pi/180 ;                 %position2 = peshawar
lat2 =  34*pi/180;                        %GC dist   = 145.8km
lon2 = (71+35/60)*pi/180;


xe1 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat1)))*cos(lat1)*cos(lon1);
ye1 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat1)))*cos(lat1)*sin(lon1);
ze1 = (r_e*(1-epsi*epsi)/sqrt(1-epsi*epsi*sin(lat1)*sin(lat1)))*sin(lat1);

xe2 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat2)))*cos(lat2)*cos(lon2);
ye2 = (r_e/sqrt(1-epsi*epsi*sin(lat1)*sin(lat2)))*cos(lat2)*sin(lon2);
ze2= (r_e*(1-epsi*epsi)/sqrt(1-epsi*epsi*sin(lat1)*sin(lat2)))*sin(lat2);

d = sqrt((xe1-xe2)*(xe1-xe2) + (ye1-ye2)*(ye1-ye2));
fprintf('total distance =  %f\n',d);


%accel1
nbr = 5/deltat;                           %changed from orignal code
lat2t = (33+49/60)*pi/180;
lon2t = (72)*pi/180;
ntax = 1000;                                  % 1.2*60/deltat                     %ntax = 720

vitax = 0; 
vftax = 60;
psi1 = 140*pi/180;




%......................................................................................................................
%AlTERATION ALGORITHM COMPUTING THE BEARING 
deg=1;
nm=60.02*deg;
m=1852*nm;

%{
%NOTE: I have added this part to compute the required heading between two points, so that the heading calculated could be used as psi1.
%NOTE: Rhumb line bearing came out to  be 86 deg...but using it as psi still didnt do any good to results. 
%NOTE: if Great Circle bearing (-86 deg) used as psi, it moved the lat lon in the right direction but very very slowly. 

%great circle computations
distancegc=2*asin(sqrt(func1(lat2t-lat1) + cos(lat1).*cos(lat2t).*func1(lon2t-lon1)));
    distancegc= distancegc*180*m/pi;
    fprintf('great circle distance=  %f\n',distancegc);
    
    bearinggc = atan2(sin(lon2t-lon1).*cos(lat2t),cos(lat1).*sin(lat2t) - sin(lat1).*cos(lat2t).*cos(lon2t-lon1))*180/pi;
    fprintf('great circle bearing=  %0.8f\n',bearinggc);
%}
    
%{
%rhumb line
     y=log(tan(lat2t./2+pi/4)./tan(lat1./2+pi/4));
  if (lat2t==lat1)
     v = cos(lat1);
   else 
     v = (lat2t-lat1)/y;
  end   
  if (lon2t > lon1)
      course = mod(180.*(atan2((lon1-lon2t),y))./pi,360);
      totaldistance= 180*m.*(sqrt(v^2*(lon2t-lon1)^2 + (lat2t-lat1)^2))./pi;
   else
      course = mod(180.*(atan2((lon1-lon2t),y))./pi,360);
      totaldistance= 180*m.*(sqrt(v.^2.*(lon1-lon2t).^2 + (lat2t-lat1).^2))./pi;
  end
     fprintf('Rhumb line distance =  %f\n',totaldistance);
     fprintf('Rhumb line course =  %0.8f\n',course);
     
    
%}

%........................................................................................................................



acceltraj(lat1,lon1,lat2t,lon2t,ntax,0,vftax,vitax,deltat,r_e,psi1,0,0,0,0,vitax*cos(psi1),vitax*sin(psi1),0,0,0,0,0,xe1,ye1,ze1);
ntotal = ntax;
trjd




nclimb = 10;                          %inorder to test the acceleration function only ,put breakpoint here
lat2cl = (33+50/60)*pi/180
lon2cl = (73)*pi/180
vicl=70;                              %QUESTION :: its never used?
vfcl = 130;
vfx = -5;
vfy = 60;
height(ntotal)=0;

deacceltraj(lata(ntotal),lona(ntotal),lat2cl,lon2cl,nclimb,ntotal,vfcl,trjd.vfin,deltat,r_e,trjd.psi,height(ntotal),1000*0.3048,-5,-3,trjd.vx,trjd.vy,trjd.theta,trjd.phi,vfx,vfy,trjd.vex,trjd.vey,trjd.vez,xea(ntotal),yea(ntotal),zea(ntotal))
ntotal = ntotal + nclimb;
trjd



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

