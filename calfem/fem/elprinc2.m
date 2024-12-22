function [magnfac]=elprinc2(ex,ey,es,plotprop,magnfac)
%elprinc2(ex,ey,es,plotpar,magnfac)
%[magnfac]=elprinc2(ex,ey,es,plotpar)
%[magnfac]=elprinc2(ex,ey,es)
%-------------------------------------------------------------
% PURPOSE 
%   Display element principal stresses as arrows for a number of  
%   2D structural elements of the same type. To display elements and 
%   nodes, use eldraw2. Supported elements are:
%
%     1) -> triangular 3 node el.    2) -> quadrilateral 4 node el. 
%
% INPUT    
%    ex,ey:.......... nen:   number of element nodes
%                     nel:   number of elements   
%    es:     element stress matrix
%
%    plotpar=[  arrowtype, arrowcolor]
%
%        arrowtype=1 -> solid       arrowcolor=1 -> white
%                  2 -> dashed                 2 -> green
%                  3 -> dotted                 3 -> yellow
%                                              4 -> red
%        
%    magnfac:  =  arrowlength / max element principal
%
%    Rem. Default is auto magnification and solid white arrows if magnfac
%         and plotpar is left out.
%         
%-------------------------------------------------------------

% LAST MODIFIED: P-E Austrell 1994-01-07 
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
 if ~((nargin==3)|(nargin==4)|(nargin==5))
    error('??? Wrong number of input arguments!')
 end
 
 a=size(ex); b=size(ey); c=size(es);
 
 if (a-b)==[0 0]
     nel=a(1);nen=a(2); 
 else
    error('??? Check size of coordinate input arguments!') 
    k
 end
 
 if ~((nen==3)|(nen==4))  
    error('Sorry, this element is currently not supported!') 
    
 end    
 
 if ~(c(1)==a(1))
    disp('??? Check size of stress input argument!')
    error('One row for each element, i.e the mean stress sig x, -y and tau xy !') 
     
 end
% 
% if ned~=nen ; 
%    disp('??? This function should be used for structure problems!')
% end
%
% ******** calculation of principal stresses and directions **********

 sigx=es(:,1); sigy=es(:,2); tau=es(:,3);
 
 ds=(sigx-sigy)/2; R=sqrt(ds.^2+tau.^2);
 
 sig1=(sigx+sigy)/2+R; sig2=(sigx+sigy)/2-R; alfa=atan2(tau,ds)/2;
 
% ************************************************************************
 dxmax=max(max(ex')-min(ex'));  dymax=max(max(ey')-min(ey'));
 
 lm=sum(sqrt(dxmax.^2+dymax.^2))/nel;
 
 sig1m=sum(sig1)/nel; sig2m=sum(sig2)/nel; sigm=max(sig1m,sig2m);

 krel=0.4;
 
 if nargin==3; 
    plotprop=[1 1];   magnfac=lm*krel/sigm;
 elseif nargin==4;
    magnfac=lm*krel/sigm;
 end

 plotpar=[plotprop 0]; s1=pltstyle(plotpar);
 
 x0=sum(ex')/nen; y0=sum(ey')/nen;
 
 la1=magnfac*sig1;  la2=magnfac*sig2;
 
 qrt=(pi/2).*ones(size(alfa));
 
% ************* plot commands *******************

 loc1=-sign(sig1);  loc2=-sign(sig2);

 fi=alfa;       arrow2(x0,y0,la1,fi,loc1,s1);
 
 fi=alfa+2*qrt;   arrow2(x0,y0,la1,fi,loc1,s1);      

 fi=alfa+qrt;   arrow2(x0,y0,la2,fi,loc2,s1);
 
 fi=alfa+3*qrt;   arrow2(x0,y0,la2,fi,loc2,s1);      
%--------------------------end--------------------------------
