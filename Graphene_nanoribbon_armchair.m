clear, clc;
format long
%%
N=8;
Numberofunit=1;
Nunit=2*N;
NTOT=Numberofunit*Nunit;
a=1.42;
%%
%%%%%Pos1 is matrix which consists of the positions (in the x-y plane)
%%%of 2N carbon atoms in a unit cell of an armchair GNR%%%%
Pos1=zeros(NTOT,2);
for m=0:Numberofunit-1
 for n=0:1
    for i=1:N
    if rem(i,2)==1
        Pos1(i+n*N+m*2*N,1)=a+n*a+m*3*a;
    end
    if rem(i,2)==0
        Pos1(i+n*N+m*2*N,1)=(a/2)+n*2*a+m*3*a;
    end 
    end
    for i=1:N
        Pos1(i+n*N+m*2*N,2)=(i-1)*(3^0.5)*a/2;
    end
 end
end
%%
%proj2, (b)
%%%%%%%%%%%%%the elements of Disij indicate distance between
%%%%%%%%%%%%%atom i of Pos1 and atom j of Pos2.%%%%%%%%%%%%
Pos2=Pos1;
Dis=zeros(NTOT,NTOT);
for i=1:NTOT
    for j=1:NTOT
        Dis(i,j)=(((Pos1(i,1)-Pos2(j,1))^2)+((Pos1(i,2)-Pos2(j,2))^2))^(0.5);
    end
end

%%
%proj2,  (c)

%%%%%%%%%%%%%the elements of  Hamij indicating the hopping parameter%%%%%%
%%%%%%%%%%%%%between atom i and atom j.%%%%%%%%%
t=-2.7;
Dis2=zeros(NTOT,NTOT);
for i=1:NTOT
    for j=1:NTOT
        if Dis(i,j)<= 2
            Dis2(i,j)=Dis(i,j);
        end
    end
end

Ham=(t/a)*Dis2;
Dis2=zeros(2*N,2*N);
for i=1:2*N
    for j=1:2*N
        if Dis(i,j)<= 2
            Dis2(i,j)=Dis(i,j);
        end
    end
end
%%
h01=zeros(2*N,2*N);
for i=1:N
    if rem(i,2)==0
        h01(i,i+N)=t;
    end
end
%%
KX=-pi/(3*a):pi/(3*a*25):pi/(3*a);
d=zeros(2*N,length(KX));
i=0;
for kx=KX
    i=i+1;
    H=Ham+h01'*exp(1j*kx*3*a)+h01*exp(-1j*kx*3*a);
    d(:,i)=eig(H);
    d=real(sort(d));  
end
%%
figure(1)
 for j=1:2*N
        plot(KX,(d(j,:)),'b','linewidth',2)
        hold on
        grid on
 end
       
  
   