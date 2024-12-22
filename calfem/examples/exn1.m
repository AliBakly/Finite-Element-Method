% example exn1 
%----------------------------------------------------------------
% PURPOSE 
%    Analysis of a plane frame using second order theory.
%----------------------------------------------------------------

% REFERENCES
%     Susanne Heyden 95-11-01 
%     Karl-Gunnar Olsson 96-01-23
%----------------------------------------------------------------
 echo off

% ----- Topology -----

 Edof=[1  1  2  3  4  5  6;
       2 10 11 12  7  8  9;
       3  4  5  6  7  8  9]; 

% ----- Element properties and global coordinates ----- 
      
 E=2.1e11;  
 A1=45.3e-4;     A2=142.8e-4;
 I1=2510e-8;	 I2=33090e-8;
 ep1=[E A1 I1];	 ep3=[E A2 I2];

 Ex=[0 0;6 6;0 6];   Ey=[0 4;0 4;4 4];

% ----- Load vector -----

 f=zeros(12,1);	
 f(4)=1000;   f(5)=-1000000;   f(8)=-1000000;

% ----- Initial values for the iteration -----

 eps=0.001;		% Error norm
 N=[0.01 0 0];		% Initial normal forces
 N0=[1 1 1];		% Normal forces of the initial former iteration
 n=0;			% Iteration counter

% ----- Iteration procedure -----

 while(abs((N(1)-N0(1))/N0(1))>eps)
   n=n+1

   K=zeros(12,12);

   Ke1=beam2g(Ex(1,:),Ey(1,:),ep1,N(1));
   Ke2=beam2g(Ex(2,:),Ey(2,:),ep1,N(2));
   Ke3=beam2g(Ex(3,:),Ey(3,:),ep3,N(3));

   K=assem(Edof(1,:),K,Ke1);
   K=assem(Edof(2,:),K,Ke2);
   K=assem(Edof(3,:),K,Ke3);

   bc=[1 0;2 0;3 0;10 0;11 0;12 0];	
   a=solveq(K,f,bc);

   Ed=extract(Edof,a);

   es1=beam2gs(Ex(1,:),Ey(1,:),ep1,Ed(1,:),N(1));
   es2=beam2gs(Ex(2,:),Ey(2,:),ep1,Ed(2,:),N(2));
   es3=beam2gs(Ex(3,:),Ey(3,:),ep3,Ed(3,:),N(3));

   N0=N
   N=[es1(1,1),es2(1,1),es3(1,1)];
   
   if(n>20)
     disp('The solution doesn''t converge')
     return
   end
 end
 echo off


