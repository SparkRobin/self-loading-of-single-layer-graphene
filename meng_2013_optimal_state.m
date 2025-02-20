% Xianhong Meng et al 2013 J. Phys. D: Appl. Phys. 46 055308
clear
clc

EI=1.4; %eV.bending stiffness
d=0.34; %nm. interlayer distance between parallel region
gamma=1.45; %binding energy. per unit area of graphene. 𝑒𝑉/𝑛𝑚^2

L=2.5; %optimal half length of curved region. nm
Ltot=11; %total length=2(L+L0).nm

%% theoretical k0,k1 for L=2.5 d=0.34

k0=1.02;
k1=1.95;
disp('With optimal state k0=1.02 k1=1.915')

f1=@(theta) 1./sqrt(k0.^2+sin(theta).*(k1.^2-k0.^2));
f2=@(theta) 1./sqrt(k0.^2-sin(theta).*(k1.^2-k0.^2));
f3=@(theta) sin(theta)./sqrt(k0.^2-sin(theta).*(k1.^2-k0.^2));
f4=@(theta) sin(theta)./sqrt(k0.^2+sin(theta).*(k1.^2-k0.^2));

Ln=integral(f1,0,pi/2)+2*integral(f2,0,asin(k0^2./(k1.^2-k0.^2)));
Ln=real(Ln)
dn=-2*integral(f3,0,asin(k0^2./(k1.^2-k0.^2)))+integral(f4,0,pi/2);
dn=dn*2;
dn=real(dn)

%Energy
f5=@(theta) sqrt(k0.^2-sin(theta).*(k1.^2-k0.^2));
f6=@(theta) sqrt(k0.^2+sin(theta).*(k1.^2-k0.^2));
Ubend=EI*(2*integral(f5,0,asin(k0^2./(k1.^2-k0.^2)))+integral(f6,0,pi/2))
Utotal=Ubend-gamma/2*(Ltot-2*Ln)

%% How to find optimal k0,k1

ii=1;
jj=1;
k0_map=0.8:0.002:1.2;
k1_map=1.3:0.002:2; % k1>k0
for k0=k0_map
    for k1=k1_map
        f1=@(theta) 1./sqrt(k0.^2+sin(theta).*(k1.^2-k0.^2));
        f2=@(theta) 1./sqrt(k0.^2-sin(theta).*(k1.^2-k0.^2));
        f3=@(theta) sin(theta)./sqrt(k0.^2-sin(theta).*(k1.^2-k0.^2));
        f4=@(theta) sin(theta)./sqrt(k0.^2+sin(theta).*(k1.^2-k0.^2));
        Ln=integral(f1,0,pi/2)+2*integral(f2,0,asin(k0^2./(k1.^2-k0.^2)));
        Ln_map(ii,jj)=real(Ln);
        dn=-2*integral(f3,0,asin(k0^2./(k1.^2-k0.^2)))+integral(f4,0,pi/2);
        dn=dn*2;
        dn_map(ii,jj)=real(dn);
        jj=jj+1;

        % if abs(Ln-L)<=0.03 && abs(dn-d)<=0.003
        %     k0
        %     k1
        %     keyboard
        % end
    end
    jj=1;
    ii=ii+1;
end

[k0_map,k1_map]=meshgrid(k1_map,k0_map);

%% Plot

surf(k0_map,k1_map,dn_map,'EdgeColor','none')
xlabel('k1')
ylabel('k0')
title('d')
hold on
z = d*ones(size(dn_map));
surf(k0_map,k1_map,z,'FaceColor','r','EdgeColor','none')
hold off

figure
surf(k0_map,k1_map,Ln_map,'EdgeColor','none')
xlabel('k1')
ylabel('k0')
title('L')
hold on
z = L*ones(size(dn_map));
surf(k0_map,k1_map,z,'FaceColor','r','EdgeColor','none')
hold off