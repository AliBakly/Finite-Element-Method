   function G=trans3(ex,ey,ez,e0)
% G=trans3(ex,ey,ez,e0)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the transformation matrix for a 3D  
%  beam element. 
% 
% INPUT:  ex = [x1 x2];     node coordinates, and 
%         ey = [y1 y2];     coordinates for a point on 
%         ez = [z1 z2];     the positive local y axis
%
%         e0 = [x0 y0 z0];  coordinate on the local y-axis
%                       
%    OUTPUT: G : transformation matrix (12 x 12)
%    
%-------------------------------------------------------------
 
% REFERENCES
%   P-E AUSTRELL 1993-09-24 
%   K Persson    1995-08-23
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology
%-------------------------------------------------------------
  
    b=[ ex(2)-ex(1); ey(2)-ey(1); ez(2)-ez(1) ];
    L=sqrt(b'*b); n1=b/L;
%
    c=[ e0(1)-ex(1); e0(2)-ey(1); e0(3)-ez(1) ];
    lc=sqrt(c'*c); n2=c/lc;
%
    n3(1)=n1(2)*n2(3)-n1(3)*n2(2);
    n3(2)=n1(3)*n2(1)-n1(1)*n2(3);
    n3(3)=n1(1)*n2(2)-n1(2)*n2(1);
%
    An=[n1';
        n2';
        n3];
%
    G=[  An     zeros(3) zeros(3) zeros(3);
       zeros(3)   An     zeros(3) zeros(3);
       zeros(3) zeros(3)   An     zeros(3);
       zeros(3) zeros(3) zeros(3)   An    ];
%
%--------------------------end--------------------------------
