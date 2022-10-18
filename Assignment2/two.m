%Code to solve the simultaneous equations(negative x direction) for using TDMA and plot the Non-Dimensional Temperature Theta vs x
clear
L=1;
Pe_c=10;
m=25;
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
disp(Theta)
X=linspace(0,-L,m+1);
plot(X,Theta)
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
xlabel("<----- X");
ylabel("Theta ---->");
title("Theta_2 vs X (Negative X Direction)");
saveas(1,"Plot-2.png")
