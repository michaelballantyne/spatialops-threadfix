clear; clc; close all;
doX = 1;
doY = 1;
doZ = 1;

test2der = 0;

if( doX )
   % x-data
   fx=mmread('fx.mm'); x=mmread('x.mm'); gx=mmread('gx.mm'); xg=mmread('xg.mm');
   df=mmread('dfdx.mm');
   subplot(3,1,1); plot(x,fx,'k.',xg,gx,'ro'); legend('orig','interp');
   subplot(3,1,2); plot(xg,cos(xg),'kx',xg,df,'rs'); legend('dx','num');
   if( test2der )
      d2f=mmread('d2fdx2.mm');
      subplot(3,1,3); plot(x,-sin(x),'kx',x,d2f,'rs'); legend('d2x','num');
   end
   title('x results');
end

if ( doY )
   % y-data
   figure;
   fy=mmread('fy.mm'); y=mmread('y.mm'); gy=mmread('gy.mm'); yg=mmread('yg.mm');
   df=mmread('dfdy.mm');
   subplot(3,1,1);  plot(y,fy,'k.',yg,gy,'ro');     legend('orig','interp');
   subplot(3,1,2);  plot(yg,cos(yg),'kx',yg,df,'rs'); legend('exact grad','num grad');
   if( test2der )
      d2f=mmread('d2fdy2.mm');
      subplot(3,1,3);  plot(y,-sin(y),'kx',y,d2f,'bs');   legend('exact div','num div');
   end
   title('y results');
end

if( doZ )
   % z-data
   figure;
   fz=mmread('fz.mm'); z=mmread('z.mm'); gz=mmread('gz.mm'); zg=mmread('zg.mm');
   df=mmread('dfdz.mm');
   subplot(3,1,1);  plot(z,fz,'k.',zg,gz,'ro');     legend('orig','interp');
   subplot(3,1,2);  plot(zg,cos(zg),'kx',zg,df,'rs'); legend('exact grad','num grad');
   if( test2der )
      d2f=mmread('d2fdz2.mm');
      subplot(3,1,3);  plot(z,-sin(z),'kx',z,d2f,'bs');   legend('exact div','num div');
   end
   title('z results');
end