%% Tight bonding method calculation: band structure of bilayer graphene
clear, clc;
format long

a1=3.28373;
a2=4.63208;
t0=1.0;  %% hopping parameter (eV)
t1=-1.62;
e=0.38*t0;
l=-0.1*t0;%% SOC strength
% Reciprocal lattice vectors
%a1 =  [1.0, 0.0];
%a2 =  [0.0, 1.0];
%r1 = [1.0,  0.59215];
%r2 = [0.5,  0.90785];
%r3=  [1.0,  0.40785];
% define k-line
Kx= pi/(a1); 
Ky= pi/(a2);
u=1;
d=-1;
%%
%Highsymmetry points generation
points=5000; % Mesh points

% Define the path from -X  to Γ to X to Y 
path = [-Kx, 0; 0, 0; Kx, 0; 0, Ky];

% % Generate additional points between -X  to Γ
 X1x_to_gammax =  linspace(path(1, 1), path(2, 1), points);
 X1y_to_gammay = linspace(path(1, 2), path(2, 2), points);

% % Generate additional points between Γ to X
 gammax_to_nXx = linspace(path(2, 1), path(3, 1), points);
 gammay_to_nXy = linspace(path(2, 2), path(3, 2), points);

% % Generate additional points between X to Y 
 nXx_to_Yx = linspace(path(3, 1), path(4, 1), points);
 nXy_to_Yy = linspace(path(3, 2), path(4, 2), points);

% % Concatenate all points to create the extended path
 kx_path = [X1x_to_gammax gammax_to_nXx nXx_to_Yx];
 ky_path = [X1y_to_gammay gammay_to_nXy nXy_to_Yy];


% % x axis for ploting band structure
 kx_plot=[linspace(0,2,length(X1x_to_gammax)) ...
          linspace(2,4,length(gammax_to_nXx))...
          linspace(4,4+(2*sqrt(2)),length(nXx_to_Yx))];
 % % x axis for ploting band structure
 ky_plot=[linspace(0,2,length(X1y_to_gammay)) ...
          linspace(2,4,length(gammay_to_nXy))...
          linspace(4,4+(2*sqrt(2)),length(nXy_to_Yy))];
%%
for count=1:1:length(kx_path)

    kxn=kx_path(count);
    kyn=ky_path(count);

    % Construct hamiltonian

    f1= 2*t0*(exp(1i*1.46234.*kyn)).*cos(kxn*(1.64186));

    f2= t1*(exp(-1i*(kyn)*(-0.8537)));
    
    f3= 2*l*sin(a1.*kxn); %%SOC function

    % Hamiltonian1=[e+(u*f3), f1,   0,   f2',        0,  0,  0,  0;
    %              f1',- e-(u*f3), f2,   0,          0,  0,  0,  0;           
    %              0,     f2',   u*f3,  f1,          0,  0,  0,  0;
    %              f2,    0,     f1',-(u*f3),        0,  0,  0,  0;  
    %              0,  0,  0,  0,                    0,  0,  0,  0;
    %              0,  0,  0,  0,                    0,  0,  0,  0;
    %              0,  0,  0,  0,                    0,  0,  0,  0;
    %              0,  0,  0,  0,                    0,  0,  0,  0
    %                                               ];

    Hamiltonian1=[e+(u*f3), f1,   0,   f2';
                 f1',- e-(u*f3), f2,   0;           
                 0,     f2',   u*f3,  f1;
                 f2,    0,     f1',-(u*f3)
           
                                                  ];

% Hamiltonian2=[   0,  0,  0,  0,              0,  0,  0,  0;
%                  0,  0,  0,  0,              0,  0,  0,  0;           
%                  0,  0,  0,  0,              0,  0,  0,  0;
%                  0,  0,  0,  0,              0,  0,  0,  0;  
%                  0,  0,  0,  0,            e+(d*f3), f1,    0, f2';
%                  0,  0,  0,  0,            f1',  -e-(d*f3), f2, 0;
%                  0,  0,  0,  0,            0,      f2',    d*f3, f1;
%                  0,  0,  0,  0,            f2,     0,    f1', -(d*f3)
%                                                   ];
Hamiltonian2=[
              e+(d*f3), f1,    0, f2';
              f1',  -e-(d*f3), f2, 0;
              0,      f2',    d*f3, f1;
              f2,     0,    f1', -(d*f3)
                                                  ];

    % get eigen values of Hamiltonian
    Eig1(count,:)=sort(real(eig(Hamiltonian1)),'ascend');
    Eig2(count,:)=sort(real(eig(Hamiltonian2)),'ascend');

end

%% Plot band structure

figure(4);
plot(kx_plot,Eig1,'color','r','linewidth',2)
hold on 
plot(kx_plot,Eig2,'color','b','linewidth',2)

set(gca,'fontsize',28)
set(gca,'XTick',[])
ylabel('E (eV)','FontSize',28)
set(gcf,'Position',[400 100 800 600])
xlim([0,max(kx_plot)])
grid on 
grid minor
y_l=1.2;
ylim([-y_l,y_l])

% Gamma points
x=2.*ones(1,points+1);
y=-y_l:2*y_l./points:y_l;
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 

% -K points
x=4.*ones(1,points+1);
y=-y_l:2*y_l./points:y_l;
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 

% Y points
x=(4+(2*sqrt(2))).*ones(1,points+1);
y=-y_l:2*y_l./points:y_l;
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 
