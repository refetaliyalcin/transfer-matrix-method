clc
clear all

lamda=(300:10:25000)'*10^-9;%meter
h_1=300*10^-9; %meter
h_2=10*10^-9; %meter
h_3=300*10^-9; %meter
ang_rad=(0:1:89)'*pi/180;

d=[NaN,h_1,h_2,h_3,NaN];

n_air=ones(length(lamda),1);
n_silica=sio2_n(lamda);
k_silica=sio2_k(lamda);
n_silver=silver_n(lamda);
k_silver=silver_k(lamda);

n=[n_air, n_silica + k_silica*1j, n_silver + k_silver*1j, n_silica + k_silica*1j, n_air];


if size(n,1) == length(lamda) 
    if size(n,2)==length(d)
        R=zeros(length(lamda),length(ang_rad));
        T=zeros(length(lamda),length(ang_rad));
        A=zeros(length(lamda),length(ang_rad));
        for j=1:length(ang_rad)
            for i=1:length(lamda)

                [R_TM,T_TM,A_TM]=jreftran_rt(lamda(i),d,n(i,:),ang_rad(j),1); %TM
                [R_TE,T_TE,A_TE]=jreftran_rt(lamda(i),d,n(i,:),ang_rad(j),0); %TE

                R(i,j) = 0.5*(R_TE + R_TM); %unpolarized
                T(i,j) = 0.5*(T_TE + T_TM); %unpolarized
                A(i,j) = 0.5*(A_TE + A_TM); %unpolarized
            end
        end
        figure
        plot(lamda*10^6,A(:,1),lamda*10^6,T(:,1),'LineWidth',2)
        xlabel('Wavelength [\mum]')
        ylabel('Spectral Absorptance and Transmittance')
        ylim([0 1])
        legend('A','T')
            
    else
        error('Array sizes are inconsistent. Thickness vector is not consistent with refractive index matrix')
    end
else
    error('Array sizes are inconsistent. Wavelength vector is not consistent with refractive index matrix ')
end