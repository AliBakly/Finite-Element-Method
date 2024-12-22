function [N,dNr,blk]=xshape(nen,gp)
% [N]=xshape(nen,gp) 
% [N,dNr,blk]=xshape(nen,gp) 
%-------------------------------------------------------------
% PURPOSE
%  Calculate the values of the interpolation functions
%  an its derivatives. Store pointers for blocks in blk.
%-------------------------------------------------------------

% REFERENCES
%  M Ristinmaa 1993-08-28
%  H Carlsson  1994-02-06
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology
%-------------------------------------------------------------
[ngp,nsd]=size(gp);

xsi=gp(:,1);
%
% 3-node element
%
 if nen==3
   N(:,1)= .5*xsi.*(xsi-1);
   N(:,2)=  1-xsi.*xsi;
   N(:,3)= .5*xsi.*(xsi+1);
   if nargout>1
     dNr(:,1)= xsi-0.5;
     dNr(:,2)= -2*xsi;
     dNr(:,3)= xsi+0.5;
   end
 end
if nsd==2
 eta=gp(:,2);
%
% 4-node element
%
 if nen==4
   N(:,1)=(1-xsi).*(1-eta);
   N(:,2)=(1+xsi).*(1-eta);
   N(:,3)=(1+xsi).*(1+eta);
   N(:,4)=(1-xsi).*(1+eta);
   N=N*.25;
   if nargout>1
     r2=ngp*nsd;
     dNr(1:2:r2,1)=-(1-eta);
     dNr(1:2:r2,2)= (1-eta);
     dNr(1:2:r2,3)= (1+eta);
     dNr(1:2:r2,4)=-(1+eta);
     dNr(2:2:r2+1,1)=-(1-xsi);
     dNr(2:2:r2+1,2)=-(1+xsi);
     dNr(2:2:r2+1,3)= (1+xsi);
     dNr(2:2:r2+1,4)= (1-xsi);
     dNr=dNr*.25;
   end
 elseif nen==8
%
% 8-node element
%
  N(:,1)=-((1-xsi).*(1-eta).*(1+xsi+eta))*.25;
  N(:,2)=-((1+xsi).*(1-eta).*(1-xsi+eta))*.25;
  N(:,3)=-((1+xsi).*(1+eta).*(1-xsi-eta))*.25;
  N(:,4)=-((1-xsi).*(1+eta).*(1+xsi-eta))*.25;
  N(:,5)=((1-xsi.*xsi).*(1-eta))*.5;
  N(:,6)=((1+xsi).*(1-eta.*eta))*.5;
  N(:,7)=((1-xsi.*xsi).*(1+eta))*.5;
  N(:,8)=((1-xsi).*(1-eta.*eta))*.5;
  if nargout>1
    r2=ngp*nsd;
    dNr(1:2:r2,1)=-(-(1-eta).*(1+xsi+eta)+(1-xsi).*(1-eta))*.25;
    dNr(1:2:r2,2)=-( (1-eta).*(1-xsi+eta)-(1+xsi).*(1-eta))*.25;
    dNr(1:2:r2,3)=-( (1+eta).*(1-xsi-eta)-(1+xsi).*(1+eta))*.25;
    dNr(1:2:r2,4)=-(-(1+eta).*(1+xsi-eta)+(1-xsi).*(1+eta))*.25;
    dNr(1:2:r2,5)=-xsi.*(1-eta);
    dNr(1:2:r2,6)=(1-eta.*eta)*.5;
    dNr(1:2:r2,7)=-xsi.*(1+eta);
    dNr(1:2:r2,8)=-(1-eta.*eta)*.5;
    dNr(2:2:r2+1,1)=-(-(1-xsi).*(1+xsi+eta)+(1-xsi).*(1-eta))*.25;
    dNr(2:2:r2+1,2)=-(-(1+xsi).*(1-xsi+eta)+(1+xsi).*(1-eta))*.25;
    dNr(2:2:r2+1,3)=-( (1+xsi).*(1-xsi-eta)-(1+xsi).*(1+eta))*.25;
    dNr(2:2:r2+1,4)=-( (1-xsi).*(1+xsi-eta)-(1-xsi).*(1+eta))*.25;
    dNr(2:2:r2+1,5)=-(1-xsi.*xsi)*.5;
    dNr(2:2:r2+1,6)=-eta.*(1+xsi);
    dNr(2:2:r2+1,7)=(1-xsi.*xsi)*.5;
    dNr(2:2:r2+1,8)=-eta.*(1-xsi);
  end
 elseif  nen==9
%
% 9-node element
%
   disp('unknown element');
   fmess;
   return
 end;
%
% solids
%
elseif nsd==3
 eta=gp(:,2);
 zet=gp(:,3);
 if nen==8
  N(:,1)=(1-xsi).*(1-eta).*(1-zet);
  N(:,2)=(1+xsi).*(1-eta).*(1-zet);
  N(:,3)=(1+xsi).*(1+eta).*(1-zet);
  N(:,4)=(1-xsi).*(1+eta).*(1-zet);
  N(:,5)=(1-xsi).*(1-eta).*(1+zet);
  N(:,6)=(1+xsi).*(1-eta).*(1+zet);
  N(:,7)=(1+xsi).*(1+eta).*(1+zet);
  N(:,8)=(1-xsi).*(1+eta).*(1+zet);
  N=N/8.;
  if nargout>1
    r2=ngp*nsd;
    dNr(1:3:r2,1)=-(1-eta).*(1-zet);
    dNr(1:3:r2,2)= (1-eta).*(1-zet);
    dNr(1:3:r2,3)= (1+eta).*(1-zet);
    dNr(1:3:r2,4)=-(1+eta).*(1-zet);
    dNr(1:3:r2,5)=-(1-eta).*(1+zet);
    dNr(1:3:r2,6)= (1-eta).*(1+zet);
    dNr(1:3:r2,7)= (1+eta).*(1+zet);
    dNr(1:3:r2,8)=-(1+eta).*(1+zet);
    dNr(2:3:r2+1,1)=-(1-xsi).*(1-zet);
    dNr(2:3:r2+1,2)=-(1+xsi).*(1-zet);
    dNr(2:3:r2+1,3)= (1+xsi).*(1-zet);
    dNr(2:3:r2+1,4)= (1-xsi).*(1-zet);
    dNr(2:3:r2+1,5)=-(1-xsi).*(1+zet);
    dNr(2:3:r2+1,6)=-(1+xsi).*(1+zet);
    dNr(2:3:r2+1,7)= (1+xsi).*(1+zet);
    dNr(2:3:r2+1,8)= (1-xsi).*(1+zet);
    dNr(3:3:r2+2,1)=-(1-xsi).*(1-eta);
    dNr(3:3:r2+2,2)=-(1+xsi).*(1-eta);
    dNr(3:3:r2+2,3)=-(1+xsi).*(1+eta);
    dNr(3:3:r2+2,4)=-(1-xsi).*(1+eta);
    dNr(3:3:r2+2,5)= (1-xsi).*(1-eta);
    dNr(3:3:r2+2,6)= (1+xsi).*(1-eta);
    dNr(3:3:r2+2,7)= (1+xsi).*(1+eta);
    dNr(3:3:r2+2,8)= (1-xsi).*(1+eta);
    dNr=dNr/8.;
  end
 end
end

if nargout>1
  % Compute pointers to block of dNr for each 
  % gauss point. 

  j=1:ngp;
  blk0=nsd*(j-1);
  for i=1:nsd
    blk(i,:)=blk0+i; 
  end
end
%--------------------------end--------------------------------
