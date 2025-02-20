% Xianhong Meng et al 2013 J. Phys. D: Appl. Phys. 46 055308
clear
clc

EI=1.4; %eV.bending stiffness
d=0.34; %nm. interlayer distance between parallel region
gamma=1.45; %binding energy. per unit area of graphene. 𝑒𝑉/𝑛𝑚^2

L=2.5; %optimal half length of curved region. nm
Ltot=11; %total length=2(L+L0).nm
k0=1.02;
k1=1.95;
%% Plot profile from eqs.(14)

disp('With optimal state k0=1.02 k1=1.95')

f1=@(theta) cos(theta)./sqrt(k0.^2-sin(theta).*(k1.^2-k0.^2));
f2=@(theta) sin(theta)./sqrt(k0.^2-sin(theta).*(k1.^2-k0.^2));
theta1=asin(k0^2./(k1.^2-k0.^2));

profile1=[];
ii=1;
for ktheta=linspace(0,theta1,1000)
xAB=real(integral(f1,0,ktheta));
yAB=d/2+real(integral(f2,0,ktheta));
profile1(ii,:)=[xAB yAB];
ii=ii+1;
end

profile2=[];
ii=1;
for ktheta=linspace(theta1,0,1000)
xBC=real(integral(f1,0,theta1)+...
    integral(f1,ktheta,theta1));
yBC=d/2+real(integral(f2,0,theta1)+...
    integral(f2,ktheta,theta1));
profile2(ii,:)=[xBC yBC];
ii=ii+1;
end

profile3=[];
ii=1;
for ktheta=linspace(0,-pi/2,1000)
xCD=real(2*integral(f1,0,theta1)-...
    integral(f1,0,ktheta));
yCD=d/2+real(2*integral(f2,0,theta1)-...
    integral(f2,0,ktheta));
profile3(ii,:)=[xCD yCD];
ii=ii+1;
end


scatter(profile1(:,1),profile1(:,2),'filled')
hold on
scatter(profile2(:,1),profile2(:,2),'filled')
scatter(profile3(:,1),profile3(:,2),'filled')
axis equal
legend('AB','BC','CD')
xlabel('x')
ylabel('y')
title('Large deformation theoretical model')
ax=gca;
ax.FontSize=15;
ax.FontName='Arial';
ax.FontWeight='bold';