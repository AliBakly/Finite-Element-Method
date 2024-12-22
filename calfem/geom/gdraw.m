function gdraw(typ,h,nr)
% gdraw(typ,h,nr)
%------------------------------------------------------------------------------
% PURPOSE
%  Draw part of the geometry
%
% INPUT:  typ  type of geometry element to be drawn ('point','line' or 'surf')
%
%         h    if h=0 the figure is erased before drawn
%
%         nr   number of the paet to be drawn, if omited all numbers
%-------------------------------------------------------------

% REFERENCES
%   M Ristinmaa 1994-05-24
%   A Magnusson 1994-05-24
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology

global CALFEM_POINTS;
global CALFEM_LINES;
global CALFEM_SURFS;
if nargin==1 h=1; end
if h==0
   clf
end
hold on
typs=typ(1:3);
if typs=='poi'
   if nargin==3
     r=find(CALFEM_POINTS(:,1)==nr);
     if isempty(r)
       disp('POINT ERROR: pointnumber not found')
     else
       plot(CALFEM_POINTS(r,2),CALFEM_POINTS(r,3),'g*');
       text(CALFEM_POINTS(r,2),CALFEM_POINTS(r,3),[' P',num2str(nr)]);
     end
   else
     [r,c]=size(CALFEM_POINTS);
     for ii=1:r
       if CALFEM_POINTS(ii,1)~=0
          gdraw('poi',1,CALFEM_POINTS(ii,1))
       end
     end
   end
elseif typs=='lin'
   if nargin==3
     r=find(CALFEM_LINES(:,1)==nr);
     if isempty(r)
       disp('LINE ERROR: linenumber not found')
       nr
     else
       p1=CALFEM_LINES(r,3); p2=CALFEM_LINES(r,4); p3=CALFEM_LINES(r,5);
         rp=[find(CALFEM_POINTS(:,1)==p1) find(CALFEM_POINTS(:,1)==p2)];
         if p3==0
           xy=[CALFEM_POINTS(rp,2) CALFEM_POINTS(rp,3)]';
           plot(xy(1,:),xy(2,:),'b-')
           xc=sum(xy')/2;
           text(xc(1),xc(2),[' L',num2str(nr)]);
          else
           rpc=[find(CALFEM_POINTS(:,1)==p3)];
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
           div=abs(CALFEM_LINES(r,2));
           df=fiTot/div;
           for ii=1:div
             fi1=df*(ii-1);
             fi2=df*(ii);
             r1=ra+(rb-ra)*fi1/fiTot;
             r2=ra+(rb-ra)*fi2/fiTot;
             plot([xc+r1*cos(fi0-fi1),xc+r2*cos(fi0-fi2)],[yc+r1*sin(fi0-fi1),yc+r2*sin(fi0-fi2)],'b-')
             if ii==round(div/2)
               text(xc+r1*cos(fi0-fi1),yc+r1*sin(fi0-fi1),[' L',num2str(nr)]);
             end
           end
         end
     end
   else 
     [r,c]=size(CALFEM_LINES);
     for i=1:r
       if CALFEM_LINES(i,2)==0
         nl=find(CALFEM_LINES(i,:)==0);
         if isempty(nl) nl=c+1; end
         for ii=3:nl-1
           gdraw('lin',1,CALFEM_LINES(i,ii))
         end
       else
         gdraw('lin',1,CALFEM_LINES(i,1))
       end
     end
   end
elseif typs=='sur'
   if nargin==3
     r=find(CALFEM_SURFS(:,1)==nr);
     if isempty(r)
       disp('SURFACE ERROR: surfacenumber not found')
     else
       c=find(CALFEM_SURFS(r,:)==0);
       if isempty(c) [rr,c]=size(CALFEM_SURFS); end
       for ii=2:c
          gdraw('lin',1,CALFEM_SURFS(r,ii))
       end
     end
   else 
     [r,c]=size(CALFEM_SURFS);
     for i=1:r
       gdraw('sur',1,CALFEM_SURFS(i,1))
     end
   end
end
%axis('square')
%axis('image')

