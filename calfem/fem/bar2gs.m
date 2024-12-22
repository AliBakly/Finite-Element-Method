function [es]=bar2gs(ex,ey,ep,ed,N)
% es=bar2gs(ex,ey,ep,ed,N)
%-------------------------------------------------------------
% PURPOSE
%  Compute the normal force in the two dimensional 
%  bar element.
%
% INPUT:  ex = [x1 x2]
%         ey = [y1 y2]         element coordinates
%
%         ep = [E A]           E : Young's modulus
%                              A : Cross section area
%
%         ed : [u1, ..., u4;   element displacement vector
%               ...........]   one row for each element
%
%         N                    normal force
%
%  OUTPUT: es = [N1 ;          normal force in element 
%                N2 ;
%               ...]           one row for each element
%-------------------------------------------------------------

% REFERENCES
%   P-E Austrell 1993-08-26
%   K Persson    1994-01-18
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology
%-------------------------------------------------------------
   E=ep(1);  A=ep(2);
%
   b=[ ey(2)-ey(1); ey(2)-ey(1) ];
   L=sqrt(b'*b);
   n=b/L;
%
   G=[n(1) n(2)  0   0   ;
     -n(2) n(1)  0   0   ;
       0    0   n(1) n(2);
       0    0  -n(2) n(1)];
%
   Kle=E*A/L*[1  0 -1  0;
              0  0  0  0;
             -1  0  1  0;
              0  0  0  0]+N/L*[0  0  0  0;
                               0  1  0 -1;
                               0  0  0  0;
                               0 -1  0  1];
%                               
   u=ed';
   P1=(E*A/L*[-1 0 1 0]*G*u)';
   es=P1;
%--------------------------end--------------------------------
