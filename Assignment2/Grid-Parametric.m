%Code to perform the grid Sensitivity Study and Parametric Study
clear
L=1;
Pe_c=10;
m=1;
Theta_M=[];
Theta_P=[];
M=[];
P=[];
p=0.1;
j=1;
%Grid Sensitivity Study
while(m<=50)
  Theta=zeros(1,m+1);
  Theta(1)=1;
  Theta(m+1)=0;
  D=(-(2+Pe_c))*ones(1,m);
  D(1)=1;
  A=(1+Pe_c)*ones(1,m);
  A(m)=0;
  A(1)=0;
  B=ones(1,m);
  B(1:2)=0;
  C=zeros(1,m);
  C(1)=1;
  C(2)=-1;
  for i=2:m
    w=B(i)/D(i-1);
    D(i)=D(i)-w*A(i-1);
    C(i)=C(i)-w*C(i-1);
  endfor
  Theta(m)=C(m)/D(m);
  for i=m-1:-1:1
    Theta(i)=(C(i)-A(i)*Theta(i+1))/D(i);
  endfor
  Theta_M(m)= Theta(2);
  M(m)=m;
  m+=1;
endwhile
%Parametric Study for cell Peclet Number (Pe_c)
while(p<=10)
  m=10;
  Pe_c=p;
  Theta=zeros(1,m+1);
  Theta(1)=1;
  Theta(m+1)=0;
  D=(-(2+Pe_c))*ones(1,m);
  D(1)=1;
  A=(1+Pe_c)*ones(1,m);
  A(m)=0;
  A(1)=0;
  B=ones(1,m);
  B(1:2)=0;
  C=zeros(1,m);
  C(1)=1;
  C(2)=-1;
  for i=2:m
    w=B(i)/D(i-1);
    D(i)=D(i)-w*A(i-1);
    C(i)=C(i)-w*C(i-1);
  endfor
  Theta(m)=C(m)/D(m);
  for i=m-1:-1:1
    Theta(i)=(C(i)-A(i)*Theta(i+1))/D(i);
  endfor
  Theta_P(j)= Theta(2);
  P(j)=p;
  p+=0.1;
  j+=1;
endwhile
%Plotting Theta vs number of nodes
figure(1)
plot(M,Theta_M)
xlabel("Number of elements (m)")
ylabel("Theta_2")
title("Grid-Sensitivity Study")
hold off
%Plotting Theta vs cell peclet number
figure(2)
plot(P,Theta_P)
xlabel("Peclet Number")
ylabel("Theta_2")
title("Parametric Study (Theta vs Pe_c)")
saveas(1,"Grid_Sensitivity.png")
saveas(2,"Parmetric.png")
