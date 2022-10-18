clear
L=1;
W=1;
m=1;
a=L/W;
A=a^2;
Bi=2;%h*L/k;
lambda=1.9;
tol=1e-4;
M=[];
Theta_L=[];
while(m<=25)
  n=m/a;
  X=linspace(0,1,m+1);
  Y=linspace(0,1,n+1);
  x=linspace(0,L,m+1);
  y=linspace(0,W,n+1);
  DX=a*1/m;
  DY=1/n;
  Theta=ones(m+1,n+1);
  maxErr=100;
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
  endwhile
  Theta_L(m)=Theta(m+1,n+1);
  M(m)=m;
  printf("m=%f\titerations=%f\tTheta=%f\n",m,iter,Theta(m+1,n+1))
  m+=1;
endwhile
plot(M,Theta_L)
xlabel("Number of nodes")
ylabel("Top-Right Tip Temperature(\Theta_m_+_1_,_n_+_1)")
title("Grid Sensitivity Study")
saveas(1,"Grid-Sensitivity.png");