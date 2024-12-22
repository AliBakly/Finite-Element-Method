function [X,Y,nodes,xcn,ycn]=nextcurved(x1,y1,x2,y2,xcn,ycn,r1,r2,rel1,rel2,ny,div,k,dofn,X,Y,nodes)

global CALFEM_POINTS
global CALFEM_LINES
global CALFEM_COORD
global CALFEM_DOF

nextpoint=size(CALFEM_POINTS,1)+1;
nextline=size(CALFEM_LINES,1)+1;
nextnode=size(CALFEM_COORD,1)+1;

d=sqrt((x2-x1)^2+(y2-y1)^2);
rp=[x2-x1 y2-y1];
rp=rp/norm(rp);
ro=[rp(2) -rp(1)];
b=(r1^2-r2^2)/2/d+d/2;
a=sqrt(abs(r1^2-b^2));
xc1=x1+b*rp(1)+a*ro(1);
yc1=y1+b*rp(2)+a*ro(2);
xc2=x1+b*rp(1)-a*ro(1);
yc2=y1+b*rp(2)-a*ro(2);
if norm([xc2-xcn;yc2-ycn])<norm([xc1-xcn;yc1-ycn])
   xcn=xc2; ycn=yc2;
else
   xcn=xc1; ycn=yc1;
end
gpoints(nextpoint,x1,y1,0,0)
gpoints(nextpoint+1,x2,y2,0,0)
gpoints(nextpoint+2,xcn,ycn,0,0)
rel=(rel2-rel1)/ny*(k-1)+rel1;
glines(nextline,div,nextpoint,nextpoint+1,nextpoint+2,rel,0)
[X(k,:), Y(k,:)]=divline(nextline);
nodes(k,2:div)=nextnode:nextnode+div-2;
CALFEM_COORD=[CALFEM_COORD; X(k,2:div)' Y(k,2:div)' zeros(div-1,1)];
newdofs=[];
for j=1:dofn
   newdofs=[newdofs; (nextnode-1)*dofn+j:dofn:(nextnode+div-2)*dofn];
end
CALFEM_DOF=[CALFEM_DOF; newdofs'];
CALFEM_POINTS(nextpoint:nextpoint+2,:)=[];
CALFEM_LINES(nextline,:)=[];
