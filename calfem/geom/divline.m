function [x,y]=divline(line)
% calculates x and y coordinates for the nodes on the
% divided line li

% REFERENCES
%   M Ristinmaa 1994-05-24
%   A Magnusson 1994-05-24
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology

global CALFEM_POINTS
global CALFEM_LINES
li=CALFEM_LINES(line,:);
if li(2)~=0    % not a combined line
   rel=li(6);
   div=li(2);
   rp=[find(CALFEM_POINTS(:,1)==li(3)); find(CALFEM_POINTS(:,1)==li(4))];
   if li(5)==0    % straigt line
      Lx=CALFEM_POINTS(rp(2),2)-CALFEM_POINTS(rp(1),2);
      Ly=CALFEM_POINTS(rp(2),3)-CALFEM_POINTS(rp(1),3);
      L=sqrt(Lx^2+Ly^2);
      if abs(rel)<=1 rel=0; end
      if rel==0
         lvect=L/div*ones(1,div);
      else
         lvect=1:(abs(rel)-1)/(div-1):abs(rel);
         lvect=L/sum(lvect)*lvect;
      end
      xs=CALFEM_POINTS(rp(1),2);
      ys=CALFEM_POINTS(rp(1),3);
      xe=CALFEM_POINTS(rp(2),2);
      ye=CALFEM_POINTS(rp(2),3);
      if sign(rel)<0
         lvect=fliplr(lvect);
      end
      x=cumsum(lvect)*Lx/L+xs;
      y=cumsum(lvect)*Ly/L+ys;
      x=[xs x];
      y=[ys y];
      x(div+1)=xe;
      y(div+1)=ye;
   else
      rel=-rel;
      rpc=find(CALFEM_POINTS(:,1)==li(5));
      x1=CALFEM_POINTS(rp(1),2); y1= CALFEM_POINTS(rp(1),3);
      x2=CALFEM_POINTS(rp(2),2); y2= CALFEM_POINTS(rp(2),3);
      xc=CALFEM_POINTS(rpc,2)  ; yc= CALFEM_POINTS(rpc,3);
      ra=sqrt((x1-xc)^2+(y1-yc)^2);
      rb=sqrt((x2-xc)^2+(y2-yc)^2);
      fi0=atan2((y1-yc),(x1-xc));
      fi1=atan2((y2-yc),(x2-xc));
      fiTot=fi0-fi1;
      if fiTot<-pi
         fiTot=fiTot+2*pi;
      elseif fiTot>pi
         fiTot=fiTot-2*pi;
      end
      trn=0;
      if rel<0 trn=1; end
      if abs(rel)<=1 
         rel=0;
      end
      if rel==0
         df(1)=0;
         df=[df fiTot/div*ones(1,div-1)];
         df=cumsum(df);
      else
         if rel<0 rel=1/rel; end
         c=2*fiTot*(1-abs(rel))/(div*(div-1)*(abs(rel)+1));
         dfi0=abs(rel)*c*(div-1)/(1-abs(rel));
         df(1)=0;
         for ii=2:div
            df(ii)=dfi0+c*(ii-2);
         end
         df=cumsum(df);
      end
      xs=x1; ys=y1; xe=x2; ye=y2;
      x=xs;
      y=ys;
      for ii=1:div-1
         fi1=df(ii);
         fi2=df(ii+1);
         r1=ra+(rb-ra)*fi1/fiTot;
         r2=ra+(rb-ra)*fi2/fiTot;
         x=[x xc+r2*cos(fi0-fi2)];
         y=[y yc+r2*sin(fi0-fi2)];
      end
      x(div+1)=xe;
      y(div+1)=ye;
   end   
else
   k=3;
   x=[]; y=[];
   while li(k) 
      [x1,y1]=divline(CALFEM_LINES(li(k),:));
      if k>3 x1(1)=[]; y1(1)=[]; end
      x=[x x1];
      y=[y y1];
      k=k+1;
    end
end
