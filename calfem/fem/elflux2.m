function [magnfac]=elflux(ex,ey,ef,plotprop,magnfac)
%elflux2(ex,ey,ef,plotprop,magnfac)
%[magnfac]=elflux2(ex,ey,ef,plotprop)
%[magnfac]=elflux2(ex,ey,ef)
%-------------------------------------------------------------
% PURPOSE 
%   Display element flow arrows for a number of 2D scalar elements of 
%   the same type. To display elements and nodes, use eldraw2. 
%   Supported elements are:
%
%     1) -> triangular 3 node el.    2) -> quadrilateral 4 node el. 
%
% INPUT    
%    ex,ey:.......... nen:   number of element nodes
%                     nel:   number of elements   
%    ef:     element flow matrix
%
%    plotprop=[  arrowtype, arrowcolor]
%
%        arrowtype=1 -> solid       arrowcolor=1 -> white
%                  2 -> dashed                 2 -> green
%                  3 -> dotted                 3 -> yellow
%                                              4 -> red
%        
%    magnfac: magnification factor  =  arrowlength / element flow
%
%    Rem. Default is auto magnification and solid white arrows if magnfac
%         and plotpar is left out.
%         
%-------------------------------------------------------------

% LAST MODIFIED: P-E Austrell 1994-01-06 
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
 if ~((nargin==3)|(nargin==4)|(nargin==5))
    error('??? Wrong number of input arguments!')
    
 end
 
 a=size(ex); b=size(ey); c=size(ef);
 
 if (a-b)==[0 0]
     nel=a(1);nen=a(2); 
 else
    error('??? Check size of coordinate input arguments!') 
    
 end

 if ~(c(1)==a(1))
    disp('??? Check size of flow input argument!')
    error('One row for each element, i.e the mean flow in x- and y-directions !') 
     
 end
 
 ned=c(2); 
 
% if ned~=nen ; 
%    disp('??? This function should be used for scalar problems!')
% end

 dxmax=max(max(ex')-min(ex')); 
 dymax=max(max(ey')-min(ey'));
 lm=sum(sqrt(dxmax.^2+dymax.^2))/nel;
 
 qm=sum(sqrt(sum((ef').^2)))/nel;

 krel=0.8;
 
 if nargin==3; 
    plotprop=[1 1];
    magnfac=lm*krel/qm;
 elseif nargin==4;
    magnfac=lm*krel/qm;
 end
% *****************************************************************
 if ~((nen==3)|(nen==4))  
     error('Sorry, this element is currently not supported!') 
      
 else
% ************* plot commands ******************* 
    plotpar=[plotprop 0];
    s1=pltstyle(plotpar);
 
    q=sqrt(sum((ef').^2));
    la=magnfac*q; fi=atan2(ef(:,2),ef(:,1)); loc=zeros(size(la));
  
    x0=sum(ex')/nen; y0=sum(ey')/nen;
 
    arrow2(x0,y0,la,fi,loc,s1)
 end        
%--------------------------end--------------------------------
