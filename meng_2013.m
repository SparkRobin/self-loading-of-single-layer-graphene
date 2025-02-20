% Xianhong Meng et al 2013 J. Phys. D: Appl. Phys. 46 055308
clear
clc

EI=1.4; %eV.bending stiffness
d=0.34; %nm. interlayer distance between parallel region
gamma=1.45; %binding energy. per unit area of graphene. 𝑒𝑉/𝑛𝑚^2
%5,6.5,8,11,12.5
Ltot=5; %total length=2(L+L0).nm
diary(sprintf('Ltot_%d.txt',Ltot));
disp('Computation begin')

%% L,d with k0,k1 scan
tic
i=1;
paras=[];
for L=1:0.1:Ltot/2
    for k0=0.001:0.002:6
        for k1=k0+0.1:0.002:6
            f1=@(theta) 1./sqrt(k0.^2+sin(theta).*(k1.^2-k0.^2));
            f2=@(theta) 1./sqrt(k0.^2-sin(theta).*(k1.^2-k0.^2));
            f3=@(theta) sin(theta)./sqrt(k0.^2-sin(theta).*(k1.^2-k0.^2));
            f4=@(theta) sin(theta)./sqrt(k0.^2+sin(theta).*(k1.^2-k0.^2));
            Ln=integral(f1,0,pi/2)+2*integral(f2,0,asin(k0^2./(k1.^2-k0.^2)));
            dn=-2*integral(f3,0,asin(k0^2./(k1.^2-k0.^2)))+integral(f4,0,pi/2);
            dn=dn*2;
            %Energy
            f5=@(theta) sqrt(k0.^2-sin(theta).*(k1.^2-k0.^2));
            f6=@(theta) sqrt(k0.^2+sin(theta).*(k1.^2-k0.^2));
            Ubend=EI*(2*integral(f5,0,asin(k0^2./(k1.^2-k0.^2)))+integral(f6,0,pi/2));
            Utotal=Ubend-gamma/2*(Ltot-2*Ln);

            if abs(Ln-L)<=0.05 && abs(dn-d)<=0.005
                break
            end
        end
        if abs(Ln-L)<=0.05 && abs(dn-d)<=0.005
            Utotal_n(i)=Utotal;
            paras(i,:)=[k0 k1];
            Utotal
            toc
            break
        end
    end
    L
    i=i+1;
end

disp('Computation Completed')

plot(1:0.1:Ltot/2,Utotal_n)
xlabel('L(nm)')
ylabel('U_{total} (eV/nm)')
legend(sprintf('L_{total}=%d mn',Ltot))
folder='data';
% Create folder (if doesn't exist)
if ~isfolder(folder)
    mkdir(folder);
end
filename = sprintf('Ltot%d.mat',Ltot);
save(fullfile(folder,filename));

filename2 = sprintf('Ltot%d.png',Ltot);
saveas(gcf,fullfile(folder,filename2),'png');

toc
disp('Data Saved Completed')
