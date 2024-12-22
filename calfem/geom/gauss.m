function  [gp,wp]=gauss(ir,nsd)
% [gp,wp]=gauss(ir,nsd)
%-------------------------------------------------------------
% PURPOSE
%  Setup vectors used for intergration of isoparametric
%  elements.
%-------------------------------------------------------------

% REFERENCES
%    M Ristinmaa 1993-08-24
%    H Carlsson  1994-02-03
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology
%-------------------------------------------------------------
% error handle
%
if nargin~=2
  error('Wrong number of input arguments ! Use: gauss(ir,nsd)');
  return
end
if ir>3
  error('Error ! used number of integration points not implemented');
  return
end
%-------------------------------------------------------------
if ir==1
  g1=0.;
  w1=2.;
elseif ir==2
  g1= 0.577350269189626;
  w1= 1.;
elseif ir==3
  g1= 0.774596699241483;
  g2= 0.;
  w1= 0.555555555555555;
  w2= 0.888888888888888;
end

if nsd==2
  gp=zeros(ir*ir,nsd);
  wp=zeros(ir*ir,nsd);
  if ir==1
    gp(1,1)=g1;
    gp(1,2)=g1;
    w(1,1)=w1;
    w(1,2)=w1;
  elseif ir==2
    gp(:,1)=[-1; 1;-1; 1]*g1;
    gp(:,2)=[-1;-1; 1; 1]*g1;
    w(:,1)=[ 1; 1; 1; 1]*w1;
    w(:,2)=[ 1; 1; 1; 1]*w1;
  elseif ir==3
    gp(:,1)=[-1; 0; 1;-1; 0; 1;-1; 0; 1]*g1;
    gp(:,1)=[ 0;-1; 0; 0; 1; 0; 0; 1; 0]*g2+gp(:,1);
    gp(:,2)=[-1;-1;-1; 0; 0; 0; 1; 1; 1]*g1;
    gp(:,2)=[ 0; 0; 0; 1; 1; 1; 0; 0; 0]*g2+gp(:,2);
    w(:,1)=[ 1; 0; 1; 1; 0; 1; 1; 0; 1]*w1;
    w(:,1)=[ 0; 1; 0; 0; 1; 0; 0; 1; 0]*w2+w(:,1);
    w(:,2)=[ 1; 1; 1; 0; 0; 0; 1; 1; 1]*w1;
    w(:,2)=[ 0; 0; 0; 1; 1; 1; 0; 0; 0]*w2+w(:,2);
  else
    disp('Used number of integration points not implemented');
    return
  end

  wp=w(:,1).*w(:,2);   

elseif nsd==3
  gp=zeros(ir*ir*ir,nsd);
  w=zeros(ir*ir*ir,nsd);
  if ir==2
    gp(:,1)=[-1; 1; 1;-1;-1; 1; 1;-1]*g1;
    gp(:,2)=[-1;-1; 1; 1;-1;-1; 1; 1]*g1;
    gp(:,3)=[-1;-1;-1;-1; 1; 1; 1; 1]*g1;
    w(:,1)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    w(:,2)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    w(:,3)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
  elseif ir==3
    I1=[-1; 0; 1;-1; 0; 1;-1; 0; 1]';
    I2=[ 0;-1; 0; 0; 1; 0; 0; 1; 0]';
    gp(:,1)=[I1 I1 I1]'*g1;
    gp(:,1)=[I2 I2 I2]'*g2+gp(:,1);
    I1=abs(I1);
    I2=abs(I2);
    w(:,1)=[I1 I1 I1]'*w1;
    w(:,1)=[I2 I2 I2]'*w2+w(:,1);
    I1=[-1;-1;-1; 0; 0; 0; 1; 1; 1]';
    I2=[ 0; 0; 0; 1; 1; 1; 0; 0; 0]';
    gp(:,2)=[I1 I1 I1]'*g1;
    gp(:,2)=[I2 I2 I2]'*g2+gp(:,2);
    I1=abs(I1);
    I2=abs(I2);
    w(:,2)=[I1 I1 I1]'*w1;
    w(:,2)=[I2 I2 I2]'*w2+w(:,2);
    I1=[-1;-1;-1;-1;-1;-1;-1;-1;-1]';
    I2=[ 0; 0; 0; 0; 0; 0; 0; 0; 0]';
    I3=abs(I1);
    gp(:,3)=[I1 I2 I3]'*g1;
    gp(:,3)=[I2 I3 I2]'*g2+gp(:,3);
    w(:,3)=[I3 I2 I3]'*w1;
    w(:,3)=[I2 I3 I2]'*w2+w(:,3);
  else
    disp('Used number of integration points not implemented');
    return
  end

  wp=w(:,1).*w(:,2).*w(:,3);   

end



%--------------------------end--------------------------------
