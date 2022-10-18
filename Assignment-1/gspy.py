import numpy as np
import matplotlib.pyplot as plt
L=1
W=1
q=50
k=2
h=1
T_f=25
m=5
n=5
x=np.linspace(0,L,m+1)
y=np.linspace(0,W,n+1)
a=1#L/W
A=a**2
X=np.linspace(0,1,m+1)
Y=np.linspace(0,1,n+1)
DX=1/m
DY=1/n
Bi=0.5#h*L/k
Theta=np.ones((m+1,n+1))
#print(Theta)
maxErr=100
tol=1e-4
iter=0
while (maxErr>tol):
    maxErr=0
    for i in range(0,m+1):
        for j in range(n+1):
            if ((i+1)==1 and (j+1)==1):
                d=((2*(DY**2)*Theta[i+1,[j]])+(2*A*(DX**2)*Theta[i,[j+1]])+((DX*DY)**2))/(2*((DY**2)+A*(DX**2)))
            elif (((i+1)==1)and((j+1)>=2 and (j+1)<=(n))):
                d=((2*(DY**2)*Theta[i+1,[j]])+(A*(DX**2)*(Theta[i,[j-1]]+Theta[i,[j+1]]))+((DX*DY)**2))/(2*((DY**2)+(A*(DX**2))))
            elif (((i+1)>=2 and (i+1)<=(m))and((j+1)==1)):
                d=(((DY**2)*(Theta[i-1,[j]]+Theta[i+1,[j]]))+(2*A*(DX**2)*Theta[i,[j+1]])+((DX*DY)**2))/(2*((DY**2)+(A*(DX**2))))
            elif (((i+1)>=2 and (i+1)<=(m))and((j+1)>=2 and (j+1)<=(n))):
                d=(((DY**2)*(Theta[i-1,[j]]+Theta[i+1,[j]]))+(A*(DX**2)*(Theta[i,[j-1]]+Theta[i,[j+1]]))+((DX*DY)**2))/(2*((DY**2)+(A*(DX**2))))
            elif (((i+1)==(m+1))and((j+1)>=2 and (j+1)<=(n))):
                d=((2*(DY**2)*Theta[i-1,[j]])+(A*(DX**2)*(Theta[i,[j-1]]+Theta[i,[j+1]]))+((DX*DY)**2))/(2*((DY**2)+(A*(DX**2))+(Bi*(DX)*(DY**2))))
            elif (((i+1)>=2 and (i+1)<=(m))and((j+1)==(n+1))):
                d=(((DY**2)*(Theta[i-1,[j]]+Theta[i+1,[j]]))+(2*A*(DX**2)*Theta[i,[j-1]])+((DX*DY)**2))/(2*((DY**2)+(A*(DX**2))+(Bi*a*(DX**2)*(DY))))
            elif (((i+1)==(m+1))and((j+1)==(n+1))):
                d=((2*(DY**2)*Theta[i-1,[j]])+(2*A*(DX**2)*Theta[i,[j-1]])+((DX*DY)**2))/(2*((DY**2)+(A*(DX**2))+(Bi*(a*(DX**2)*(DY)+(DX)*(DY**2)))))
            elif (((i+1)==1)and((j+1)==(n+1))):
                d=((2*(DY**2)*Theta[i+1,[j]])+(2*A*(DX**2)*Theta[i,[j-1]])+((DX*DY)**2))/(2*((DY**2)+(A*(DX**2))+(Bi*a*(DX**2)*(DY))))
            elif (((i+1)==(m+1))and((j+1)==1)):
                d=((2*(DY**2)*Theta[i-1,[j]])+(2*A*(DX**2)*Theta[i,[j+1]])+((DX*DY)**2))/(2*((DY**2)+(A*(DX**2))+(Bi*(DX)*(DY**2))))
            e=((abs(d-Theta[i,[j]]))/((d+Theta[i,[j]])/2))*100
            if e>maxErr:
                maxErr=e
            Theta[i,[j]]=d
    iter+=1
print("Converged!!\nNumber of Iterations: {}".format(iter))
print(Theta)
T = T_f+((q*L/k)*Theta)
print(T)
plt.contourf(X,Y,Theta)
plt.show()
