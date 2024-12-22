% example exs4 
%----------------------------------------------------------------
% PURPOSE 
%    Analysis of a plane truss.
%----------------------------------------------------------------

% REFERENCES
%     P-E Austrell 1994-03-08 
%     K-G Olsson 1995-09-28
%----------------------------------------------------------------
 echo on 

%----- Topology matrix Edof -------------------------------------

 Edof=[1   1  2  5  6;
       2   3  4  7  8;
       3   5  6  9 10;
       4   7  8 11 12;
       5   7  8  5  6;
       6  11 12  9 10;
       7   3  4  5  6;
       8   7  8  9 10;
       9   1  2  7  8;
      10   5  6 11 12];
 
%----- Stiffness matrix K and load vector f ---------------------

 K=zeros(12); 
 f=zeros(12,1);	f(11)=0.5e6*sin(pi/6);	f(12)=-0.5e6*cos(pi/6);

%----- Element properties ---------------------------------------

 A=25.0e-4;	E=2.1e11;	ep=[E A];	

%----- Element coordinates --------------------------------------

 ex1=[0 2];   ex2=[0 2];   ex3=[2 4];   ex4=[2 4];   ex5=[2 2];
 ex6=[4 4];   ex7=[0 2];   ex8=[2 4];   ex9=[0 2];   ex10=[2 4];
 
 ey1=[2 2];   ey2=[0 0];   ey3=[2 2];   ey4=[0 0];   ey5=[0 2];
 ey6=[0 2];   ey7=[0 2];   ey8=[0 2];   ey9=[2 0];   ey10=[2 0];
 
%----- Element stiffness matrices  ------------------------------

 Ke1=bar2e(ex1,ey1,ep);	 Ke2=bar2e(ex2,ey2,ep);
 Ke3=bar2e(ex3,ey3,ep);	 Ke4=bar2e(ex4,ey4,ep);
 Ke5=bar2e(ex5,ey5,ep);	 Ke6=bar2e(ex6,ey6,ep);
 Ke7=bar2e(ex7,ey7,ep);	 Ke8=bar2e(ex8,ey8,ep);
 Ke9=bar2e(ex9,ey9,ep);	 Ke10=bar2e(ex10,ey10,ep);
 
%----- Assemble Ke into K ---------------------------------------

 K=assem(Edof(1,:),K,Ke1);	K=assem(Edof(2,:),K,Ke2); 
 K=assem(Edof(3,:),K,Ke3);	K=assem(Edof(4,:),K,Ke4);
 K=assem(Edof(5,:),K,Ke5);	K=assem(Edof(6,:),K,Ke6);
 K=assem(Edof(7,:),K,Ke7);	K=assem(Edof(8,:),K,Ke8); 
 K=assem(Edof(9,:),K,Ke9);	K=assem(Edof(10,:),K,Ke10);
 
%----- Solve the system of equations ----------------------------

 bc= [1 0;2 0;3 0;4 0];   
 a=solveq(K,f,bc)

%----- Element forces -------------------------------------------

 ed1=extract(Edof(1,:),a);	ed2=extract(Edof(2,:),a);
 ed3=extract(Edof(3,:),a);	ed4=extract(Edof(4,:),a);
 ed5=extract(Edof(5,:),a);	ed6=extract(Edof(6,:),a);
 ed7=extract(Edof(7,:),a);	ed8=extract(Edof(8,:),a);
 ed9=extract(Edof(9,:),a);	ed10=extract(Edof(10,:),a);


 N1=bar2s(ex1,ey1,ep,ed1)
 N2=bar2s(ex2,ey2,ep,ed2)
 N3=bar2s(ex3,ey3,ep,ed3)
 N4=bar2s(ex4,ey4,ep,ed4)
 N5=bar2s(ex5,ey5,ep,ed5)
 N6=bar2s(ex6,ey6,ep,ed6)
 N7=bar2s(ex7,ey7,ep,ed7)
 N8=bar2s(ex8,ey8,ep,ed8)
 N9=bar2s(ex9,ey9,ep,ed9)
 N10=bar2s(ex10,ey10,ep,ed10)

%---------------------------- end -------------------------------
 echo off
