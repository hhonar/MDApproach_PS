%% Simulation 6:  Comparison at various noise levels
% Purpose:  Looking at the performance of the various MD-based-PS at
% different noise levels
% Please note that the dependencies need to be used appropriately
% Written by Hamed Honari @ 2021-22 ALL RIGHTS RESERVED


clear;clc;
% Characteristics 
TR = 2;                                                 % Repetition time
fs = 1/TR;                                              % Sampling frequency
t = 0:1/fs:668-1/fs;
% Defining the freq. component of the signal
f = 0.05; % freq. component of the signal x1 & x2
N = 1000;  % number of repetition (realizations)
winLen = 30;   % window sizes for the Windowed Phase Sync. Measures


smltn = input('Which simulation do you want to run: 2. Ramp 3. Sigmoid?');
switch smltn
    case 2
        delphi = 4*pi/334.*(t-334).*(t-334>=0);
    case 3
        delphi = 2*pi./(1+exp(-0.01*(t-334)));
end
x = cos(2*pi*f*t);                                          % first signal
y = cos(2*pi*f*t + delphi);                                % second signal

noisevar = [1 4 10]; % various noise levels

%% BEMD
for q = 1:length(noisevar)
    for m = 1:N
        noise = mvnrnd([0 0],[noisevar(q) 0;0 noisevar(q)],length(t))';
        ex = noise(1,:);
        ey = noise(2,:);
        XN = x + ex;
        YN = y + ey;
        % Approach 1: create a complex signal, then decompose it 
        Z = XN + 1i*YN;                                            % complex signal
        imf{m} = bemd(Z);                                          % finding the IMFs
        % Finding the mean frequency of real and imag part of the complex imf of
        % complex signal
        mfreqX{m} = meanfreq(real(imf{m}'),fs);                                        
        mfreqY{m} = meanfreq(imag(imf{m}'),fs);

        % Phase Synchronization analysis:
        [~,ind_IMFX(m)] = min(abs(mfreqX{m} - f));
        [~,ind_IMFY(m)] = min(abs(mfreqY{m} - f));
        Data{m} = [real(imf{m}(ind_IMFX(m),:));imag(imf{m}(ind_IMFY(m),:))];
        H = hilbert(Data{m}');
        sigphase = angle(H);
        DELPHI{m} = sigphase(:,1)-sigphase(:,2);
        CCORSW{m,q} = circularslidingwindow(sigphase(:,1),sigphase(:,2),winLen,'vonmises');
        COSDELPHI1{m,q} = cos(DELPHI{m});
    end
end


% Display the Phase Sync Measures
figure(1);
for q=1:length(noisevar)
    subplot(6,3,q);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(cat(2,CCORSW{:,q}),2),0.95.*std(cat(2,CCORSW{:,q}),1,2), '-r','alpha','nan','remove');box on;
    outlinebounds(hl1,hp1)
    set([hl1],'LineWidth',2)
    legend([hp1],strcat(['Window length  = ' num2str(winLen)]),'Location','southeast');
    xlabel('time [s]','interpreter','latex');
    ylabel('CIRC','interpreter','latex');
    title('(a)','interpreter','latex');xlim([0 450]);ylim([-1 1])
    subplot(6,3,q+length(noisevar));hold on;[hl1 hp1]=boundedline((1:length(t)),mean(cat(2,COSDELPHI1{:,q}),2),0.95.*std(cat(2,COSDELPHI1{:,q}),1,2)+eps, '-m','alpha','nan','remove');box on;
    outlinebounds(hl1,hp1)
    set([hl1],'LineWidth',2)
    legend('$cos(\Delta\Phi[t])$','Location','Best');
    xlabel('time [s]','interpreter','latex');
    ylabel('CRP','interpreter','latex');
    title('(b)','interpreter','latex');xlim([0 450]);ylim([-1 1])
end


%% na-MEMD
clear COSDELPHI1 CCORSW
Data{1} = [x;y]';
indx = nchoosek(1:size(Data{1},2),2);
for q=1:length(noisevar)
    for m =1:N
        noise = mvnrnd([0 0],[noisevar(q) 0;0 noisevar(q)],length(t))';
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
        ndir = 4*n_channel_na; 
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
            CCORSW{m,q}(:,i) = circularslidingwindow(sigphase(:,1),sigphase(:,2),winLen,'vonmises');
            COSDELPHI1{m,q}(:,i) = cos(DELPHI{m}); 
        end
    end
end    

% Display the Phase Sync Measures
figure(1);
for q=1:length(noisevar)
    subplot(6,3,2*length(noisevar) + q);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(cat(2,CCORSW{:,q}),2),0.95.*std(cat(2,CCORSW{:,q}),1,2), '-r','alpha','nan','remove');box on;
    outlinebounds(hl1,hp1)
    set([hl1],'LineWidth',2)
    legend([hp1],strcat(['Window length  = ' num2str(winLen)]),'Location','southeast');
    xlabel('time [s]','interpreter','latex');
    ylabel('CIRC','interpreter','latex');
    title('(a)','interpreter','latex');xlim([0 450]);ylim([-1 1])
    subplot(6,3,3*length(noisevar) + q);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(cat(2,COSDELPHI1{:,q}),2),0.95.*std(cat(2,COSDELPHI1{:,q}),1,2)+eps, '-m','alpha','nan','remove');box on;
    outlinebounds(hl1,hp1)
    set([hl1],'LineWidth',2)
    legend('$cos(\Delta\Phi[t])$','Location','Best');
    xlabel('time [s]','interpreter','latex');
    ylabel('CRP','interpreter','latex');
    title('(b)','interpreter','latex');xlim([0 450]);ylim([-1 1])
end

%% MVMD
clear COSDELPHI1 CCORSW mfreq DELPHI data
K = 2; % num of modes apriori from na-MEMD - there could be other options to decide o

for q = 1:length(noisevar)
    for m =1:N
        noise = mvnrnd([0 0],[noisevar(q) 0;0 noisevar(q)],length(t))';
        ex = noise(1,:);
        ey = noise(2,:);
        XN = x + ex;
        YN = y + ey;
        % Tuning parameters - note that \tau = 0 regulates the Lagrangian
        % Multipliers in the case of a noisy signal (see the our paper - simulation 6 and discussion)
        % parameter for the noise 1: [u, u_hat, omega] = MVMD([XN;YN], 1000, 0, K, 1, 0, 1e-15);
        %                         4: [u, u_hat, omega] = MVMD([XN;YN], 5000, 0, K, 1, 0, 1e-15);
        %                        10: [u, u_hat, omega] = MVMD([XN;YN], 4500, 0, K, 1, 0, 1e-15);
        % Using multivariate variational mode decomposition
        [u, u_hat, omega] = MVMD([XN;YN], 1000, 0, K, 1, 0, 1e-20);
         % Finding the mean frequency of each IMFs
        for i=1:K
                mfreq{m}(i) = meanfreq(u(i,:,1),fs);                                        
        end
        [~,ind] = min(abs(mfreq{m}' - f));
        for i = 1:1:size(indx,1)
            dat{m} = [u(ind(indx(i,1)),:,1);u(ind(indx(i,1)),:,2)];
            H = hilbert(dat{m}');
            sigphase = angle(H);
            DELPHI{m} = sigphase(:,1)-sigphase(:,2);
            CCORSW{m,q} = circularslidingwindow(sigphase(:,1),sigphase(:,2),winLen,'vonmises');%,'option','window','winsize',w(1));
            COSDELPHI1{m,q} = cos(DELPHI{m});
        end
    end
end


% Display the Phase Sync Measures
figure(1);
for q=1:length(noisevar)
    subplot(6,3,4*length(noisevar) + q);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(cat(2,CCORSW{:,q}),2),0.95.*std(cat(2,CCORSW{:,q}),1,2), '-r','alpha','nan','remove');box on;
    outlinebounds(hl1,hp1)
    set([hl1],'LineWidth',2)
    legend([hp1],strcat(['Window length  = ' num2str(winLen)]),'Location','southeast');
    xlabel('time [s]','interpreter','latex');
    ylabel('CIRC','interpreter','latex');
    title('(a)','interpreter','latex');xlim([0 450]);ylim([-1 1])
    subplot(6,3,5*length(noisevar) + q);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(cat(2,COSDELPHI1{:,q}),2),0.95.*std(cat(2,COSDELPHI1{:,q}),1,2)+eps, '-m','alpha','nan','remove');box on;
    outlinebounds(hl1,hp1)
    set([hl1],'LineWidth',2)
    legend('$cos(\Delta\Phi[t])$','Location','Best');
    xlabel('time [s]','interpreter','latex');
    ylabel('CRP','interpreter','latex');
    title('(b)','interpreter','latex');xlim([0 450]);ylim([-1 1])
end



