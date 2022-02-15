%% Simulation 4:  The multicomponent signals (x: monocomponent, y: multicomponent)
% Purpose:  The purpose of theis code is to look at the mode mixing and
% ability to untangle the frequencies closely placed in signal y in the PS
% analysis in various measures.
% Please note that the dependencies need to be used appropriately
% Written by Hamed Honari @ 2021-22

clear;clc;
% Characteristics 
TR = 2;                                                 % Repetition time
fs = 1/TR;                                              % Sampling frequency
t = 0:1/fs:668-1/fs;
% Defining the freq. component of the signal
f = 0.05; % freq. component of the signal x1
N = 1000;  % number of repetition of the simulation (realizations)
w = [30 60 120];   % window sizes for the Windowed Phase Sync. Measures tWPS

smltn = input('Which simulation do you want to run: 2. Ramp 3. Sigmoid?')
switch smltn
    case 2
        delphi = 4*pi/334.*(t-334).*(t-334>=0);
    case 3
        delphi = 2*pi./(1+exp(-0.01*(t-334)));
end
x = cos(2*pi*f*t);                                          % first signal  (pure tone - ie monocomponent)
y = cos(2*pi*f*t + delphi) + cos(2*pi*f*1.1*t + delphi);    % second signal (multicomponent)


%% BEMD
for m = 1:N
noise = mvnrnd([0 0],[1 0;0 1],length(t))';
ex = noise(1,:);
ey = noise(2,:);
XN = x + ex;
YN = y + ey;
% Approach 1: create a complex signal, then decompose it 
Z = XN + 1i*YN;                                            % complex signal
imf{m} = bemd(Z);                                           % finding the IMFs
% Finding the mean frequency of real and imag part of the complex imf of
% complex signal
mfreqX{m} = meanfreq(real(imf{m}'),fs);                                        
mfreqY{m} = meanfreq(imag(imf{m}'),fs);
% Phase Synchronization analysis:
[~,ind_IMFX(m)] = min(abs(mfreqX{m} - f));
[~,ind_IMFY(m)] = min(abs(mfreqY{m} - f));
Data{m} = [real(imf{m}(ind_IMFX(m),:));imag(imf{m}(ind_IMFY(m),:))];
DatComplex{m} = imf{m}(ind_IMFX(m),:);  % using two IMFX and Y might for the index redundant
H = hilbert(Data{m}');
sigphase = angle(H);
DELPHI{m} = sigphase(:,1)-sigphase(:,2);
CCORSW{1}(:,m) = circularslidingwindow(sigphase(:,1),sigphase(:,2),w(1),'vonmises');%,'option','window','winsize',w(1));
CCORSW{2}(:,m) = circularslidingwindow(sigphase(:,1),sigphase(:,2),w(2),'vonmises');%,'option','window','winsize',w(2));
CCORSW{3}(:,m) = circularslidingwindow(sigphase(:,1),sigphase(:,2),w(3),'vonmises');%,'option','window','winsize',w(3));
COSDELPHI1(:,m) = cos(DELPHI{m});
end


Dat{1} = real(cat(1,DatComplex{1,1:N}));
Dat{2} = real(cat(1,DatComplex{1,1:N}));


% Display the Intrinsic Mode Functions
figure;
subplot(1,2,1);[hl1 hp1]=boundedline((1:length(t)),mean(Dat{1}),0.95.*std(Dat{1},1,1), '-r','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
set([hl1],'LineWidth',2)
xlabel('t');ylabel(strcat(['IMF_X ',num2str(1)]));title(strcat(['f_{mean} = ',num2str(meanfreq(mean(Dat{1}),fs),3), ' Hz']));
subplot(1,2,2);[hl1 hp1]=boundedline((1:length(t)),mean(Dat{2}),0.95.*std(Dat{2},1,1), '-b','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
set([hl1],'LineWidth',2)
xlabel('t');ylabel(strcat(['IMF_Y ',num2str(2)]));title(strcat(['f_{mean} = ',num2str(meanfreq(mean(Dat{2}),fs),3), ' Hz']));

% Display the Phase Sync Measures
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure;
subplot(3,1,1);hold on;plot((1:length(t)),delphi,'k','LineWidth',1.5);box on;
xlabel('time [s]','interpreter','latex');
ylabel('$\Delta\Phi[t] \quad [rad]$','interpreter','latex');
title('(a)','interpreter','latex');
legend('$\Delta\Phi[t]$','Location','northwest');xlim([0 450]);
subplot(3,1,2);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(CCORSW{1},2),0.95.*std(CCORSW{1},1,2), '-r','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
subplot(3,1,2);hold on;[hl2 hp2]=boundedline((1:length(t)),mean(CCORSW{2},2),0.95.*std(CCORSW{2},1,2), '-b','alpha','nan','remove');box on;
outlinebounds(hl2,hp2)
subplot(3,1,2);hold on;[hl3 hp3]=boundedline((1:length(t)),mean(CCORSW{3},2),0.95.*std(CCORSW{3},1,2), '-g','alpha','nan','remove');box on;
outlinebounds(hl3,hp3)
set([hl1 hl2 hl3],'LineWidth',2)
legend([hp1 hp2 hp3],strcat(['Window length  = ' num2str(w(1))]),strcat(['Window length  = ' num2str(w(2))]),strcat(['Window length  = ' num2str(w(3))]),'Location','southeast');
xlabel('time [s]','interpreter','latex');
ylabel('$\rho_{circ}$','interpreter','latex');
title('(b)','interpreter','latex');xlim([0 450]);ylim([-1 1])
subplot(3,1,3);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(COSDELPHI1,2),0.95.*std(COSDELPHI1,1,2)+eps, '-m','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
set([hl1],'LineWidth',2)
legend('$cos(\Delta\Phi[t])$','Location','Best');
xlabel('time [s]','interpreter','latex');
ylabel('$\vartheta$','interpreter','latex');
title('(c)','interpreter','latex');xlim([0 450]);ylim([-1 1])


%% na-MEMD
clear COSDELPHI1 CCORSW mfreq
Data{1} = [x;y]';


indx = nchoosek(1:size(Data{1},1),2);

w = [30 60 120];   % window sizes for the Windowed Phase Sync. Measures
for m =1:N
    noise = mvnrnd([0 0],[1 0;0 1],length(t))';
    ex = noise(1,:);
    ey = noise(2,:);
    XN = x + ex;
    YN = y + ey;
    Data{m} = [XN;YN]';
    % assigning the criterias
    stp_crit = 'stop';
    stp_vec = [0.3 0.3 0.3];
    mode = 'na_fix';
    intensity_noise = 0.75; 
    n_channel_na = size(Data{m},2);  
    ndir = 8*n_channel_na; % At least twice the number of channels (multivariate signals)
    imf = namemd(Data{m}, ndir, stp_crit, stp_vec, mode, intensity_noise, n_channel_na);
    % Finding the mean frequency of each IMFs
    for i=1:size(imf,1)
        for j = 1:size(imf{1},1)
            mfreq{m}(i,j) = meanfreq(imf{i}(j,:),fs);                                        
        end
    end
    
    
    % Phase Synchronization analysis:
    [~,ind] = min(abs(mfreq{m}' - f));
    for i = 1:1:size(indx,1)
        dat{m,i} = [imf{indx(i,1)}(ind(indx(i,1)),:);imf{indx(i,2)}(ind(indx(i,2)),:)];
        H = hilbert(dat{m,i}');
        sigphase = angle(H);
        DELPHI{m} = sigphase(:,1)-sigphase(:,2);
        CCORSW{1}(:,m) = circularslidingwindow(sigphase(:,1),sigphase(:,2),w(1),'vonmises');%,'option','window','winsize',w(1));
        CCORSW{2}(:,m) = circularslidingwindow(sigphase(:,1),sigphase(:,2),w(2),'vonmises');%,'option','window','winsize',w(2));
        CCORSW{3}(:,m) = circularslidingwindow(sigphase(:,1),sigphase(:,2),w(3),'vonmises');%,'option','window','winsize',w(3));
        COSDELPHI1(:,m) = cos(DELPHI{m}); 
    end

end
    

numModes = size(imf{1},1);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure;
subplot(3,1,1);hold on;plot((1:length(t)),delphi,'k','LineWidth',1.5);box on;
xlabel('time [s]','interpreter','latex');
ylabel('$\Delta\Phi[t] \quad [rad]$','interpreter','latex');
title('(a)','interpreter','latex');
legend('$\Delta\Phi[t]$','Location','northwest');xlim([0 450]);
subplot(3,1,2);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(CCORSW{1},2),0.95.*std(CCORSW{1},1,2), '-r','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
subplot(3,1,2);hold on;[hl2 hp2]=boundedline((1:length(t)),mean(CCORSW{2},2),0.95.*std(CCORSW{2},1,2), '-b','alpha','nan','remove');box on;
outlinebounds(hl2,hp2)
subplot(3,1,2);hold on;[hl3 hp3]=boundedline((1:length(t)),mean(CCORSW{3},2),0.95.*std(CCORSW{3},1,2), '-g','alpha','nan','remove');box on;
outlinebounds(hl3,hp3)
set([hl1 hl2 hl3],'LineWidth',2)
legend([hp1 hp2 hp3],strcat(['Window length  = ' num2str(w(1))]),strcat(['Window length  = ' num2str(w(2))]),strcat(['Window length  = ' num2str(w(3))]),'Location','southeast');
xlabel('time [s]','interpreter','latex');
ylabel('$\rho_{circ}$','interpreter','latex');
title('(b)','interpreter','latex');xlim([0 450]);ylim([-1 1])
subplot(3,1,3);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(COSDELPHI1,2),0.95.*std(COSDELPHI1,1,2)+eps, '-m','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
set([hl1],'LineWidth',2)
legend('$cos(\Delta\Phi[t])$','Location','Best');
xlabel('time [s]','interpreter','latex');
ylabel('$\vartheta$','interpreter','latex');
title('(c)','interpreter','latex');xlim([0 450]);ylim([-1 1])


%% MVMD 
clear COSDELPHI1 CCORSW mfreq Data 
Data{1} = [x;y];


indx = nchoosek(1:size(Data{1},1),2);
k = numModes;             % apriori information from the MEMD 
for m =1:N
    noise = mvnrnd([0 0],[1 0;0 1],length(t))';
    ex = noise(1,:);
    ey = noise(2,:);
    XN = x + ex;
    YN = y + ey;
    Data{m} = [XN;YN];
    % Using multivariate variational mode decomposition
    [u, u_hat, omega] = MVMD(Data{m}, 1000, 0.01, k, 1, 0, 1e-9);
    for i = 1:k
        mfreq{m}(1,i) = meanfreq(u(i,:,1),fs);
        mfreq{m}(2,i) = meanfreq(u(i,:,2),fs);
    end

     [~,ind] = min(abs(mfreq{m}'-f));
     dat{m} = [u(ind(indx(1,1)),:,1);u(ind(indx(1,2)),:,2)];
     H = hilbert(dat{m}');
     sigphase = angle(H);
     DELPHI{m} = sigphase(:,1)-sigphase(:,2);
     CCORSW{1}(:,m) = circularslidingwindow(sigphase(:,1),sigphase(:,2),w(1),'vonmises');%,'option','window','winsize',w(1));
     CCORSW{2}(:,m) = circularslidingwindow(sigphase(:,1),sigphase(:,2),w(2),'vonmises');%,'option','window','winsize',w(2));
     CCORSW{3}(:,m) = circularslidingwindow(sigphase(:,1),sigphase(:,2),w(3),'vonmises');%,'option','window','winsize',w(3));
     COSDELPHI1(:,m) = cos(DELPHI{m}); 
end
    




set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure;
subplot(3,1,1);hold on;plot((1:length(t)),delphi,'k','LineWidth',1.5);box on;
xlabel('time [s]','interpreter','latex');
ylabel('$\Delta\Phi[t] \quad [rad]$','interpreter','latex');
title('(a)','interpreter','latex');
legend('$\Delta\Phi[t]$','Location','northwest');xlim([0 450]);
subplot(3,1,2);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(CCORSW{1},2),0.95.*std(CCORSW{1},1,2), '-r','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
subplot(3,1,2);hold on;[hl2 hp2]=boundedline((1:length(t)),mean(CCORSW{2},2),0.95.*std(CCORSW{2},1,2), '-b','alpha','nan','remove');box on;
outlinebounds(hl2,hp2)
subplot(3,1,2);hold on;[hl3 hp3]=boundedline((1:length(t)),mean(CCORSW{3},2),0.95.*std(CCORSW{3},1,2), '-g','alpha','nan','remove');box on;
outlinebounds(hl3,hp3)
set([hl1 hl2 hl3],'LineWidth',2)
legend([hp1 hp2 hp3],strcat(['Window length  = ' num2str(w(1))]),strcat(['Window length  = ' num2str(w(2))]),strcat(['Window length  = ' num2str(w(3))]),'Location','southeast');
xlabel('time [s]','interpreter','latex');
ylabel('$\rho_{circ}$','interpreter','latex');
title('(b)','interpreter','latex');xlim([0 450]);ylim([-1 1])
subplot(3,1,3);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(COSDELPHI1,2),0.95.*std(COSDELPHI1,1,2)+eps, '-m','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
set([hl1],'LineWidth',2)
legend('$cos(\Delta\Phi[t])$','Location','Best');
xlabel('time [s]','interpreter','latex');
ylabel('$\vartheta$','interpreter','latex');
title('(c)','interpreter','latex');xlim([0 450]);ylim([-1 1])


