function [Ex,Ey,Edof,bc,last_x,last_y]=mesh
clf
t=0.02;
r=1;
nsec=31;
fi=40*pi/180;
dfi=fi/nsec;
ri=r-t/2;
rm=r;
ry=r+t/2;

coord(1:3,:)=[rm*cos(fi+pi/2) rm*sin(fi+pi/2)
				ri*cos(fi-dfi+pi/2) ri*sin(fi-dfi+pi/2)
				ry*cos(fi-dfi+pi/2) ry*sin(fi-dfi+pi/2)];

Edof(1,:)=[1 1 2 3 4	5 6];			
y=-1;
for i=1:nsec-1
	y=y+3;
	coord(i*3+1:i*3+3,:)=[rm*cos(fi-i*dfi*2+pi/2) rm*sin(fi-i*dfi*2+pi/2)
								 ri*cos(fi-(i*2+1)*dfi+pi/2) ri*sin(fi-(i*2+1)*dfi+pi/2)
					 		    ry*cos(fi-(i*2+1)*dfi+pi/2) ry*sin(fi-(i*2+1)*dfi+pi/2)];
	
	Edof(i*4-2:i*4+1,:)=[i*4-2 y*2-1 y*2 (y+2)*2-1 (y+2)*2 (y+1)*2-1 (y+1)*2
						      i*4-1 y*2-1 y*2 (y+3)*2-1 (y+3)*2 (y+2)*2-1 (y+2)*2                    
						      i*4 (y+1)*2-1 (y+1)*2 (y+2)*2-1 (y+2)*2 (y+4)*2-1 (y+4)*2
						      i*4+1   (y+2)*2-1 (y+2)*2 (y+3)*2-1 (y+3)*2 (y+4)*2-1 (y+4)*2];
								
	if (fi-(i*2+1)*dfi)==0;
		last_x=(y+4)*2-1
		last_y=(y+4)*2 
	end 
end	
y=y+3;
i=i+1;
coord(i*3+1,:)=[rm*cos(fi-i*dfi*2+pi/2) rm*sin(fi-i*dfi*2+pi/2)];
Edof(i*4-2,:)=[i*4-2 y*2-1 y*2 (y+2)*2-1 (y+2)*2 (y+1)*2-1 (y+1)*2];			


bc=[1 0;2 0;(y+2)*2-1 0; (y+2)*2 0] 

 
 Dof=[1:2:max(max(Edof(2:end,:)))-1;2:2:max(max(Edof(2:end,:)))]';
 
 
 [Ex,Ey]=coordxtr(Edof,coord,Dof,3);
 
 %eldraw2(Ex,Ey,[1 2 0]);	
 %eldraw2(Ex,Ey,[1 2 0],Edof(:,1));	
 %grid
