function G=trans2(ex,ey)
% G=trans2(ex,ey)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the transformation matrix for a 2D beam 
%  element. 
%
% INPUT:  ex = [x1 x2]
%         ey = [y1 y2]        node coordinates 
%                       
% OUTPUT: G : transformation matrix (6 x 6) 
%-------------------------------------------------------------
 
%    REFERENCES
%      P-E AUSTRELL 1993-08-27
%      K Persson    1995-08-23
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology
%-------------------------------------------------------------
  b=[ ex(2)-ex(1); ey(2)-ey(1) ];
  L=sqrt(b'*b);   n=b/L;

  G=[n(1) n(2)  0    0    0   0;
    -n(2) n(1)  0    0    0   0;
      0    0    1    0    0   0;
      0    0    0   n(1) n(2) 0;
      0    0    0  -n(2) n(1) 0;
      0    0    0    0    0   1];
           
%--------------------------end--------------------------------

