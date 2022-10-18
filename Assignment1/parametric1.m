%Code to perform the parametric study with respect to the Biot Number and plot the 4 corner Non-Dimensional Temperatures vs the Biot Number
clear
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
A=a^2;
DX=a*1/m;
DY=1/n;
Bi=0.1;
Theta=ones(m+1,n+1);
BI=[];
Th_1=[];
Th_2=[];
Th_3=[];
Th_4=[];
lambda=1.9;
k=1;
while (Bi<=10)
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
  endwhile
  Bi;
  printf("Bi= %f\t[%f,%f,%f,%f]\n",Bi,Theta(1,1),Theta(m+1,1),Theta(m+1,n+1),Theta(1,n+1))
  BI(k)=Bi;
  Th_1(k)=Theta(1,1);
  Th_2(k)=Theta(m+1,1);
  Th_3(k)=Theta(m+1,n+1);
  Th_4(k)=Theta(1,n+1);
  Bi+=0.1;
  k+=1;
endwhile
subplot(2,2,1)
plot(BI,Th_4,"+")
xlabel("Biot Number")
ylabel("Theta_1_,_n_+_1")
title("Top Left Corner")
subplot(2,2,2)
plot(BI,Th_3,"+")
xlabel("Biot Number")
ylabel("Theta_m_+_1_,_n_+_1")
title("Top Right Corner")
subplot(2,2,3)
plot(BI,Th_1,"+")
xlabel("Biot Number")
ylabel("Theta_1_,_1")
title("Bottom Left Corner")
subplot(2,2,4)
plot(BI,Th_2,"+")
xlabel("Biot Number")
ylabel("Theta_m_+_1_,_1")
title("Bottom Right Corner")
saveas(1,"Theta-Biot.png");
