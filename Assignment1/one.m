%Code for solving the Simultaneous Equations using G-S algorithm and output a contour plot of the Non-dimensional Theta on the x-y plane.
clear
%Input Parameters
h=100;
k=50;
q=500;
T_f=20;
L=1;
W=1;
t=1;
a=L/W;
nodes=15;
if (L<=W)
  m=nodes;
  n=m/a;
else
  n=nodes;
  m=a*n;
endif
%Non-Dimensional Parameters
A=a^2;
X=linspace(0,1,m+1);
Y=linspace(0,1,n+1);
x=linspace(0,L,m+1);
y=linspace(0,W,n+1);
DX=A*1/m;
DY=1/n;
Dx=L/m;
Dy=W/n;
Bi=h*L/k;
lambda=1.9;
Theta=ones(m+1,n+1);
maxErr=100;
tol=1e-4;
iter=0;
while(maxErr>tol)
  maxErr=0;
  for i=1:m+1
    for j=1:n+1
      if (i==1 && j==1)
          d=((2*(DY^2)*Theta(i+1,j))+(2*A*(DX^2)*Theta(i,j+1))+((DX*DY)^2))/(2*((DY^2)+A*(DX^2)));
        elseif ((i==1)&&(j>=2 && j<=(n)))
          d=((2*(DY^2)*Theta(i+1,j))+(A*(DX^2)*(Theta(i,j-1)+Theta(i,j+1)))+((DX*DY)^2))/(2*((DY^2)+(A*(DX^2))));
        elseif ((i>=2 && i<=(m))&&(j==1))
          d=(((DY^2)*(Theta(i-1,j)+Theta(i+1,j)))+(2*A*(DX^2)*Theta(i,j+1))+((DX*DY)^2))/(2*((DY^2)+(A*(DX^2))));
        elseif ((i>=2 && i<=(m))&&(j>=2 && j<=(n)))
          d=(((DY^2)*(Theta(i-1,j)+Theta(i+1,j)))+(A*(DX^2)*(Theta(i,j-1)+Theta(i,j+1)))+((DX*DY)^2))/(2*((DY^2)+(A*(DX^2))));
        elseif ((i==(m+1))&&(j>=2 && j<=(n)))
          d=((2*(DY^2)*Theta(i-1,j))+(A*(DX^2)*(Theta(i,j-1)+Theta(i,j+1)))+((DX*DY)^2))/(2*((DY^2)+(A*(DX^2))+(Bi*(DX)*(DY^2))));
        elseif ((i>=2 && i<=(m))&&(j==(n+1)))
          d=(((DY^2)*(Theta(i-1,j)+Theta(i+1,j)))+(2*A*(DX^2)*Theta(i,j-1))+((DX*DY)^2))/(2*((DY^2)+(A*(DX^2))+(Bi*a*(DX^2)*(DY))));
        elseif ((i==(m+1))&&(j==(n+1)))
          d=((2*(DY^2)*Theta(i-1,j))+(2*A*(DX^2)*Theta(i,j-1))+((DX*DY)^2))/(2*((DY^2)+(A*(DX^2))+(Bi*(a*(DX^2)*(DY)+(DX)*(DY^2)))));
        elseif ((i==1)&&(j==(n+1)))
          d=((2*(DY^2)*Theta(i+1,j))+(2*A*(DX^2)*Theta(i,j-1))+((DX*DY)^2))/(2*((DY^2)+(A*(DX^2))+(Bi*a*(DX^2)*(DY))));
        elseif ((i==(m+1))&&(j==1))
          d=((2*(DY^2)*Theta(i-1,j))+(2*A*(DX^2)*Theta(i,j+1))+((DX*DY)^2))/(2*((DY^2)+(A*(DX^2))+(Bi*(DX)*(DY^2))));
        endif
        d=lambda*d+((1-lambda)*Theta(i,j));
        e=((abs(d-Theta(i,j)))/((d+Theta(i,j))/2))*100;
        if e>maxErr
          maxErr=e;
        endif
        Theta(i,j)=d;
    endfor
  endfor
  iter+=1;
  if rem(iter,1000)==0
    printf("Iteration %d: Max Error=%d\n",iter,maxErr)
  endif
endwhile
printf("Converged!\nNumber of Iterations: %d\n",iter)
contour(x,y,Theta')
[C,H]=contour(x,y,Theta');
title("Contour plot of \Theta on X vs Y")
zlabel("Theta")
xlabel("X ----->")
ylabel("Y ----->")
clabel(C,H)
saveas(1,"Contour-Plot.png")
T=(q*(L^2)*(Theta)/k) + T_f;
%Calculation of heat convected from the walls
q_r=h*Dy*(T(m+1,1:n)-T_f);
q_t=h*Dx*(T(1:m,n+1)-T_f);
q_l=h*((Dx+Dy)/2)*(T(m+1,n+1)-T_f);
c=sum(q_r);
v=sum(q_t);
sum=c+v+q_l
Q=q*L*W
%Comparing with given value of heat generation and calculating error
error=((sum-Q)/Q)*100
