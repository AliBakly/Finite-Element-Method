% example exs5 
%----------------------------------------------------------------
% PURPOSE 
%    Analysis of a plane frame.
%----------------------------------------------------------------

% REFERENCES
%     G"oran Sandberg 94-03-08 
%     Karl-Gunnar Olsson 95-09-28
%     Anders Olsson 99-03-01
%----------------------------------------------------------------
 echo on

%----- Topology -------------------------------------------------

 Edof=[1  1  2  3  4  5  6;
       2 10 11 12  7  8  9;
       3  4  5  6  7  8  9];  
      

%----- Stiffness matrix K and load vector f ---------------------

 K=zeros(12);	f=zeros(12,1);	f(4)=1000;

%----- Element stiffness and element load matrices  -------------

 A1=45.3e-4;     A2=142.8e-4;
 I1=2510e-8;	 I2=33090e-8;
 E=2.1e11;  

 ep1=[E A1 I1];	 ep3=[E A2 I2];
 ex1=[0 0];      ex2=[6 6];	ex3=[0 6];
 ey1=[0 4];      ey2=[0 4];	ey3=[4 4];
 eq1=[0 0];
 eq2=[0 0];
 eq3=[0 -75000];

 Ke1=beam2e(ex1,ey1,ep1);
 Ke2=beam2e(ex2,ey2,ep1);
 [Ke3,fe3]=beam2e(ex3,ey3,ep3,eq3);

%----- Assemble Ke into K ---------------------------------------

 K=assem(Edof(1,:),K,Ke1);
 K=assem(Edof(2,:),K,Ke2);
 [K,f]=assem(Edof(3,:),K,Ke3,f,fe3);

%----- Solve the system of equations and compute reactions ------

 bc=[1 0;2 0;3 0;10 0;11 0;12 0];	
 a=solveq(K,f,bc);

%----- Section forces -------------------------------------------

 Ed=extract(Edof,a);

 [es1,edi1,eci1]=beam2s(ex1,ey1,ep1,Ed(1,:),eq1,20)
 [es2,edi2,eci2]=beam2s(ex2,ey2,ep1,Ed(2,:),eq2,20)
 [es3,edi3,eci3]=beam2s(ex3,ey3,ep3,Ed(3,:),eq3,20)
 
%----- Draw normal force diagram --------------------------------
 
 figure(1)
 magnfac=eldia2(ex1,ey1,es1(:,1),eci1);
 magnitude=[3e5 0.5 0];
 eldia2(ex1,ey1,es1(:,1),eci1,magnfac);
 eldia2(ex2,ey2,es2(:,1),eci2,magnfac);
 eldia2(ex3,ey3,es3(:,1),eci3,magnfac,magnitude);
 axis([-1.5 7 -0.5 5.5])
 
%----- Draw shear force diagram ---------------------------------
 
 figure(2)
 magnfac=eldia2(ex3,ey3,es3(:,2),eci3);
 magnitude=[3e5 0.5 0];
 eldia2(ex1,ey1,es1(:,2),eci1,magnfac);
 eldia2(ex2,ey2,es2(:,2),eci2,magnfac);
 eldia2(ex3,ey3,es3(:,2),eci3,magnfac,magnitude);
 axis([-1.5 7 -0.5 5.5])

%----- Draw moment diagram --------------------------------------
 
 figure(3)
 magnfac=eldia2(ex3,ey3,es3(:,3),eci3);
 magnitude=[3e5 0.5 0];
 eldia2(ex1,ey1,es1(:,3),eci1,magnfac);
 eldia2(ex2,ey2,es2(:,3),eci2,magnfac);
 eldia2(ex3,ey3,es3(:,3),eci3,magnfac,magnitude);
 axis([-1.5 7 -0.5 5.5])

%------------------------ end -----------------------------------
 echo off
