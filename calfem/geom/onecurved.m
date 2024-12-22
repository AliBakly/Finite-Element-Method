function [X,Y,nodes]=onecurved(X,Y,nodes,borders,dofn,close,oppos,left,right)

global CALFEM_POINTS
global CALFEM_LINES
% coordinates for corners and center point for closest line
x11=CALFEM_POINTS(CALFEM_LINES(borders(close),3),2);
y11=CALFEM_POINTS(CALFEM_LINES(borders(close),3),3);
x21=CALFEM_POINTS(CALFEM_LINES(borders(close),4),2);
y21=CALFEM_POINTS(CALFEM_LINES(borders(close),4),3);
xcn=CALFEM_POINTS(CALFEM_LINES(borders(close),5),2);
ycn=CALFEM_POINTS(CALFEM_LINES(borders(close),5),3);
a1=sqrt((x11-xcn)^2+(y11-ycn)^2);
b1=sqrt((x21-xcn)^2+(y21-ycn)^2);   % radius for first and second point

% coordinates for corners and center point for opposit line
x12=CALFEM_POINTS(CALFEM_LINES(borders(oppos),3),2);
y12=CALFEM_POINTS(CALFEM_LINES(borders(oppos),3),3);
x22=CALFEM_POINTS(CALFEM_LINES(borders(oppos),4),2);
y22=CALFEM_POINTS(CALFEM_LINES(borders(oppos),4),3);

l1=sqrt((x12-x11)^2+(y12-y11)^2);
l2=sqrt((x22-x21)^2+(y22-y21)^2);
div=CALFEM_LINES(borders(close),2);
rel1=CALFEM_LINES(borders(close),6);
rel2=CALFEM_LINES(borders(oppos),6);
ny=CALFEM_LINES(borders(left),2);
for k=2:ny
   x1=X(k,1);
   y1=Y(k,1);
   x2=X(k,div+1);
   y2=Y(k,div+1);
   d1=sqrt((x1-x11)^2+(y1-y11)^2);
   d2=sqrt((x2-x21)^2+(y2-y21)^2);
   inc=3;   % constant which determines how fast radius increase
   r1=a1*l1^inc/(l1-d1)^inc;
   r2=b1*l2^inc/(l2-d2)^inc;
   [X,Y,nodes]=nextcurved(x1,y1,x2,y2,xcn,ycn,r1,r2,rel1,rel2,ny,div,k,dofn,X,Y,nodes);
end
