clear; clc; close all;
doX = 1;
doY = 1;
doZ = 1;

if( doX )
   % x-data
   fx=mmread('fx.mm'); x=mmread('x.mm'); gx=mmread('gx.mm'); xg=mmread('xg.mm');
   df=mmread('dfdx.mm');
   d2f=mmread('d2fdx2.mm');
   subplot(3,1,1); plot(x,fx,'k.-',xg,gx,'ro'); legend('orig','interp');
   subplot(3,1,2); plot(xg,cos(xg),'kx',xg,df,'rs'); legend('dx','num');
   subplot(3,1,3); plot(x,-sin(x),'kx',x,d2f,'rs'); legend('d2x','num');
   title('x results');
end

if ( doY )
   % y-data
   figure;
   fy=mmread('fy.mm'); y=mmread('y.mm'); gy=mmread('gy.mm'); yg=mmread('yg.mm');
   df=mmread('dfdy.mm');
   plot(y,fy,'k.-',yg,gy,'ro',yg,cos(yg),'kx',yg,df,'rs');
   legend('orig','interp','exact grad','num grad');
   title('y results');
end

if( doZ )
   % z-data
   figure;
   fz=mmread('fz.mm'); z=mmread('z.mm'); gz=mmread('gz.mm'); zg=mmread('zg.mm');
   df=mmread('dfdz.mm');
   plot(z,fz,'k.-',zg,gz,'ro',zg,cos(zg),'kx',zg,df,'rs');
   legend('orig','interp','exact grad','num grad');
   title('z results');
end