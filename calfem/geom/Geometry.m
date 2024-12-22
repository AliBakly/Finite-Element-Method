clear all
warning off MATLAB:break_outside_of_loop
element='qu8'
nod_elem=8;

dr=2;
r=20;                                   % notch radius (mm)
w=100;                                      % hor. elem. dim.  (mm)
h=w*dr;                                     % layer width  (mm)
el2=16*2*2;
el1=40*2*2;
el1b=8*2;                                      % vert. elem. number. 

dinc=16e-4/0.004*0;

%dinc=0.1e-4/0.004;

%dinc=8e-4/0.004;

% clear old geom

gnew('geom');
clf;


% set points

pnum=[0 0   0 0 0 0 0   0]';
x=   [0 w-r w w w 0 w   0]';
y=   [0 0   r w 0 w 2*w 2*w]';
z=   [0 0   0 0 0 0 0   0]';

gpoints(pnum,x,y,z);

% draw lines

lnum=[0 0 0 0 0 0 0]';
eldiv=[el1 el1 el2 el2 el1b el2 el1b]';
p1=[1 3 4 6 4 7 8]';
p2=[2 4 6 1 7 8 6]';

l1=[0 el1 1 2
    0 el1 3 4
    0 el2 4 6
    0 el2 6 1
    0 el1 4 7
    0 el2 7 8
    0 el1 8 6];
glines(lnum,eldiv,p1,p2);

glinec(0,3,4);
glines(0,el2*2,2,3,5);

% create surfaces

gsurfs(0,1,9,2,8);
gsurfs(0,3,5,6,7);

% mesh

gmesh('surf',1,element,2);
gmesh('surf',2,element,2);

[edof,coord,dof]=getgeom;


bc=gboundary('line',1,2);     %y-riktning längst ner             %  
bc=gboundary('line',4,1,0,bc);  % xrikt-P3-P4           %  
bc=gboundary('line',7,1,0,bc); 	% xrit P4-P7
bc=gboundary('line',6,2,1,bc); 	% xrit P4-P7

bcarc=gboundary('line',5,1,0,bc);             %  

%bc=gboundary('line',6,2,dinc,bc);          %  
load=gboundary('line',6,2,1);
load2=gboundary('line',1,2,1);
bccond=bc;

%finc=gcload('point',3,2,1);

[ex,ey]=coordxtr(edof,coord,dof,nod_elem);
bcc=gboundary('line',6,2);
bcdof=bcc(:,1);

gnew('mesh')
gmesh('surf',1,element,2);
gmesh('surf',2,element,2);
[c,Ens1,d]=getgeom;
%[Edof,coord,Dof]=getgeom;

%[Ex,Ey]=coordxtr(Edof,coord,Dof,3);
%eldraw2(Ex,Ey,[1 4 1])

[Ex,Ey]=coordxtr(edof,coord,dof,nod_elem);
nelm=length(edof(:,1))


if 1==1

a(1)=0;
Node(1,:)=coord(1,1:2);
n=1;
for i=2:length(coord)
 flag=0;
    for j=1:nelm
    t=find(edof(j,2:end)==i*2);
    if length(t)>0
       flag=1;
   end
 end
   if flag==1
      n=n+1;
      Node(n,:)=coord(i,1:2);
      a(i)=a(i-1);
   else
   a(i)=a(i-1)+1;
   end   
end

element=[ones(nelm,1) edof(:,3)/2 edof(:,5)/2 edof(:,7)/2 ...
         edof(:,9)/2 edof(:,11)/2 edof(:,13)/2 edof(:,15)/2 edof(:,17)/2];
for i=1:nelm
   for j=2:9
      element(i,j)=element(i,j)-a(element(i,j));
   end
end
ndof=2*length(Node)
nbc=length(bc)
for i=1:nbc
   bc(i,1)=bc(i,1)-a(ceil(bc(i,1)/2))*2;
end
for i=1:length(load)
   loaddof(i)=load(i,1)-a(ceil(load(i,1)/2))*2;
end
nldof=length(loaddof)
for i=1:length(load2)
   loaddof2(i)=load2(i,1)-a(ceil(load2(i,1)/2))*2;
end
nldof2=length(loaddof2)
save('../geom.mat','Node','ndof','nelm','nbc','bc','element','loaddof','nldof','loaddof2','nldof2');

end
