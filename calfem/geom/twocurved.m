function [X,Y,nodes]=twocurved(X,Y,nodes,borders,dofn,close,oppos,left,right)
% make mesh for two curved lines (close and oppos)

global CALFEM_LINES
global CALFEM_POINTS
global CALFEM_COORD
global CALFEM_DOF

% coordinates for corners and center point for closest line
x11=CALFEM_POINTS(CALFEM_LINES(borders(close),3),2);
y11=CALFEM_POINTS(CALFEM_LINES(borders(close),3),3);
x21=CALFEM_POINTS(CALFEM_LINES(borders(close),4),2);
y21=CALFEM_POINTS(CALFEM_LINES(borders(close),4),3);
xc1=CALFEM_POINTS(CALFEM_LINES(borders(close),5),2);
yc1=CALFEM_POINTS(CALFEM_LINES(borders(close),5),3);
a1=sqrt((x11-xc1)^2+(y11-yc1)^2);
b1=sqrt((x21-xc1)^2+(y21-yc1)^2);   % radius for first and second point
%pause
% coordinates for corners and center point for opposite line
x12=CALFEM_POINTS(CALFEM_LINES(borders(oppos),3),2);
y12=CALFEM_POINTS(CALFEM_LINES(borders(oppos),3),3);
x22=CALFEM_POINTS(CALFEM_LINES(borders(oppos),4),2);
y22=CALFEM_POINTS(CALFEM_LINES(borders(oppos),4),3);
xc2=CALFEM_POINTS(CALFEM_LINES(borders(oppos),5),2);
yc2=CALFEM_POINTS(CALFEM_LINES(borders(oppos),5),3);
a2=sqrt((x12-xc2)^2+(y12-yc2)^2);
b2=sqrt((x22-xc2)^2+(y22-yc2)^2);
l1=sqrt((x12-x11)^2+(y12-y11)^2);
l2=sqrt((x22-x21)^2+(y22-y21)^2);
div=CALFEM_LINES(borders(close),2);
rel1=CALFEM_LINES(borders(close),6);
rel2=CALFEM_LINES(borders(oppos),6);
ny=CALFEM_LINES(borders(left),2);
nextnode=size(CALFEM_COORD,1)+1;
nextpoint=size(CALFEM_POINTS,1)+1;
nextline=size(CALFEM_LINES,1)+1;
% check that lines curved right way
v1=[x21-x11 y21-y11]';
n1=[v1(2) -v1(1)]';
c1=[xc1-x11 yc1-y11]';
if c1'*n1>0
   xcn=xc1;
   ycn=yc1;
else
   xcn=xc2;
   ycn=yc2;
end
for k=2:ny
   x1=X(k,1);
   y1=Y(k,1);
   x2=X(k,div+1);
   y2=Y(k,div+1);
   d1=sqrt((x1-x11)^2+(y1-y11)^2);
   d2=sqrt((x2-x21)^2+(y2-y21)^2);
   r1=a1+(a2-a1)*d1/l1;
   r2=b1+(b2-b1)*d2/l2;   % radius to first and second point;
   [X,Y,nodes,xcn,ycn]=nextcurved(x1,y1,x2,y2,xcn,ycn,r1,r2,rel1,rel2,ny,div,k,dofn,X,Y,nodes);
end
