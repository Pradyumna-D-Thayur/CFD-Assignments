clear
a=0.5;
nodes=15;
m=nodes;
n=nodes;
Theta=ones(m+1,n+1);
lambda=1.9;
Bi=2;
As=[];
Th_1=[];
Th_2=[];
Th_3=[];
Th_4=[];
k=1;
while (a<=10)
  A=a^2;
  DX=1/m;
  DY=1/n;
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
  printf("a= %f\t[%f,%f,%f,%f]\n",a,Theta(1,1),Theta(m+1,1),Theta(m+1,n+1),Theta(1,n+1))
  As(k)=a;
  Th_1(k)=Theta(1,1);
  Th_2(k)=Theta(m+1,1);
  Th_3(k)=Theta(m+1,n+1);
  Th_4(k)=Theta(1,n+1);
  a=k;
  k+=1;
endwhile
subplot(2,2,1)
plot(As,Th_4,"+")
xlabel("Aspect Ratio")
ylabel("Theta_1_,_n_+_1")
title("Top Left Corner")
subplot(2,2,2)
plot(As,Th_3,"+")
xlabel("Aspect Ratio")
ylabel("Theta_m_+_1_,_n_+_1")
title("Top Right Corner")
subplot(2,2,3)
plot(As,Th_1,"+")
xlabel("Aspect Ratio")
ylabel("Theta_1_,_1")
title("Bottom Left Corner")
subplot(2,2,4)
plot(As,Th_2,"+")
xlabel("Aspect Ratio")
ylabel("Theta_m_+_1_,_1")
title("Bottom Right Corner")
saveas(1,"Theta-Aspect.png");