%%% Simulation 5. MULTIVARIATE SIGNAL
% This script is written as part of the supplementary materials for the
% paper.  The purpose of this simulation is to compare various MD
% algorithms discussed in terms of the multivariate signal.

%% Case 5: multivariate signal using Bivariate EMD
% Characteristics 
TR = 2;                                                                     % Repetition time
fs = 1/TR;                                                                  % Sampling frequency
t = 0:1/fs:1000-1/fs;

% Defining the freq. component of the signal
% freq. component of the signals to create (0.05 also corresponds to the peak...
% of the PSD of HCP data if we look at the PSD of various IMFs and also in fmri...
% analysis often this corresponds to the center freq of the bandpass which is used 0.01-0.1)
f = 0.05;
phi1 = pi*((t-100>0).*(t-250<0) + (t-300>0).*(t-500<0) + (t-600>0).*(t-800<0));
phi2 = pi*((t-100>0).*(t-250<0) +                      - (t-600>0).*(t-800<0));
phi3 = pi*(                     - (t-300>0).*(t-500<0) - (t-600>0).*(t-800<0));
x1 = cos(2*pi*f*t+phi1);                                                    % first signal
x2 = cos(2*pi*f*t +phi2);                                                   % second signal
x3 = 2*cos(2*pi*f*t +phi3);                                                 % thrid signal

figure;subplot(3,1,1);plot(phi1,'k');ylim([-1.3*pi 1.3*pi]);ylabel('\phi_{x_1}')
subplot(3,1,2);plot(phi2,'r');ylim([-1.3*pi 1.3*pi]);ylabel('\phi_{x_2}')
subplot(3,1,3);plot(phi3,'b');ylim([-1.3*pi 1.3*pi]);ylabel('\phi_{x_3}')
xlabel('t')
subplot(3,1,1);title('Ground truth phases of the signals to generate multivariate signal')

%%%%%%%% ground truth states and transitions - note fs is considered 
grdidx = zeros(size(t));
grdidx(1:50) = 2;
grdidx(125:150) = 2;
grdidx(250:500) = 2;
grdidx(50:125) = 3;
grdidx(150:250) = 1;
grdstate{1}(:,:,1) = 0.98*[1 -1 1;-1 1 -1;1 -1 1];
grdstate{1}(:,:,2) = 0.98*[1 1 1;1 1 1;1 1 1];
grdstate{1}(:,:,3) = 0.98*[1 1 -1;1 1 -1;-1 -1 1];

% creating a trivariate signal using x_i's
Data{1} = [x1;x2;x3]';
winLen = 30;   % the window length for WPS approach

indx = nchoosek(1:size(Data{1},2),2);
N=1000; % number of realizations 
nS = 3; % number of states


%% BEMD 
for m =1:N
    noise = mvnrnd([0 0 0],[1 0 0;0 1 0;0 0 1],length(t))';
    ex1 = noise(1,:);
    ex2 = noise(2,:);
    ex3 = noise(3,:);
    XN1 = x1 + ex1;
    XN2 = x2 + ex2;
    XN3 = x3 + ex3;
    Data{m} = [XN1;XN2;XN3];
for i = 1:1:size(indx,1)
    Z{i} = Data{m}(indx(i,1),:) + 1i*Data{m}(indx(i,2),:);
%    Z{i} = XN + 1i*YN;                                            % complex signal
    imf{m,i} = bemd(Z{i});                                           % finding the IMFs
    
    
    % Finding the mean frequency of real and imag part of the complex imf of
    % complex signal
    mfreqX{m,i} = meanfreq(real(imf{m,i}'),fs);                                        
    mfreqY{m,i} = meanfreq(imag(imf{m,i}'),fs);

    % Phase Synchronization analysis:

    [~,ind_IMFX(m,i)] = min(abs(mfreqX{m,i} - f));
    [~,ind_IMFY(m,i)] = min(abs(mfreqY{m,i} - f));
    dat{m,i} = [real(imf{m,i}(ind_IMFX(m,i),:));imag(imf{m,i}(ind_IMFY(m,i),:))];
    H = hilbert(dat{m,i}');
    sigphase = angle(H);
    DELPHI = sigphase(:,1)-sigphase(:,2);
    CCORSW{m}(:,i) = circularslidingwindow(sigphase(:,1),sigphase(:,2),winLen,'vonmises');
    COSDELPHI1{m}(:,i) = cos(DELPHI);
end
    

end
%% ---------------------------------------------------
%% ---------------------------------------------------

% k-means clustering of the matrices of phase synch. (number of clusters
% using DBI) upperbound for DBI was set to nS 
[popidx{1},popCorr{1},~,~,~,DBI(1)] = mykmeans(cat(1,CCORSW{:}),nS,nS);
[popidx{2},popCorr{2},~,~,~,DBI(2)] = mykmeans(cat(1,COSDELPHI1{:}),nS,nS);

%%%%% for each simulation, match them with the population state --- or the
%%%%% ground truth
for i = 1:N
% k-means clustering of the matrices of phase synch.
[idx{1,i},Corr{1,i}] = mykmeans(CCORSW{i},nS,nS);
[idx{2,i},Corr{2,i}] = mykmeans(COSDELPHI1{i},nS,nS);
[~,~,mCorr{1,i},midx{1,i}] = matchstatesClutster(grdstate{1},Corr{1,i},idx{1,i},1);
[~,~,mCorr{2,i},midx{2,i}] = matchstatesClutster(grdstate{1},Corr{2,i},idx{2,i},1);
end

for i = 1:N
    for j = 1:nS 
        ccstate{j}(:,:,i) = mCorr{1,i}(:,:,j);
        crpstate{j}(:,:,i) = mCorr{2,i}(:,:,j);
    end
end

for i = 1:nS
    meanstate{i,1} = mean(ccstate{i},3);
    meanstate{i,2} = mean(crpstate{i},3);
end

count = 0;
placeholder = [1 2 3;7 8 9];
figure;
h = tight_subplot(3, nS, [.001 .001],[.01 .001],[.05 .05]);
labels = {'CIRC','CRP'};
set(h(4:6),'visible','off')
for j = 1:nS
    for i = 1:2
        count = count + 1;
        axes(h(placeholder(count)));gsplot(meanstate{j,i});
        axis square;
        set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w')
        if i == 1
         title(strcat(['State ' num2str(j)]));
        end
        if j == 1
            ylabel(labels{i},'interpreter','latex');
        end
    end
    
end
colorbar(h(8),'Position',[0.93 0.168 0.022 0.7]); caxis([-1 1])


for i = 1:N
    MIDX1(:,i) = midx{1,i};
    MIDX2(:,i) = midx{2,i};
end

for i = 1:size(MIDX1,1)
    flag1(i) = sum( MIDX1(i,:) == grdidx(i))/N;
    flag2(i) = sum( MIDX2(i,:) == grdidx(i))/N;
end

for j = 1:nS
for i = 1:size(MIDX1,1)
    flags1(i,j) = sum( MIDX1(i,:) == j)/N;
    flags2(i,j) = sum( MIDX2(i,:) == j)/N;
end
end




figure;subplot(3,1,1);plot(grdidx,'r','LineWidth',2);ylim([0 4]);xlabel('t');
ylabel('State #');title('Ground truth state')
subplot(3,1,2);plot(flag1,'g','LineWidth',2);
title('Accuracy of correctly classifying the state');ylabel('Accuracy');xlabel('t')
subplot(3,1,3);bar(flags1,'stacked')
title('Stacked bar of the classification accuracy across time');
legend('State 1','State 2','State 3');
ylabel('State Accuracy Proportion');xlabel('t')




figure;imagesc(MIDX1');colormap winter;
imAlpha=ones(size(MIDX1'));
imAlpha(isnan(MIDX1'))=0;
imagesc(MIDX1','AlphaData',imAlpha);
set(gca,'color',1*[1 1 1]);
title('State transitions - BEMD-based-PS using tWPS: CIRC');
xlabel('time [s]');ylabel('realizations')
cb = colorbar();
% Set color labels (one for each row in RGB)
label = 1:nS; 
caxis([1,numel(label)])
cb.YTick = 1 : nS;
labelChar = label;
cb.TickLabels = labelChar(1:end-1);
cb.FontSize = 12; 

figure;imagesc(MIDX2');colormap winter
title('State transitions - BEMD-based-PS using IPS: CRP');
xlabel('time [s]');ylabel('realizations');
cb = colorbar();
% Set color labels (one for each row in RGB)
label = 1:nS; 
caxis([1,numel(label)])
cb.YTick = 1 : nS;
labelChar = label;
cb.TickLabels = labelChar(1:end-1);
cb.FontSize = 12; 


%% na-MEMD
clear DELPHI
indx = nchoosek(1:size(Data{1},2),2);
N=1000;nS = 3;

for m =1:N
    noise = mvnrnd([0 0 0],[1 0 0;0 1 0;0 0 1],length(t))';
    ex1 = noise(1,:);
    ex2 = noise(2,:);
    ex3 = noise(3,:);
    XN1 = x1 + ex1;
    XN2 = x2 + ex2;
    XN3 = x3 + ex3;
    Data{m} = [XN1;XN2;XN3]';
    % assigning the criterias
    stp_crit = 'stop';
    stp_vec = [0.3 0.3 0.3];
    mode = 'na_fix';
    intensity_noise = 0.75; 
    n_channel_na = size(Data{m},2);  
    ndir = 2*n_channel_na; 
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
        CCORSW{m}(:,i) = circularslidingwindow(sigphase(:,1),sigphase(:,2),winLen,'vonmises');
        COSDELPHI1{m}(:,i) = cos(DELPHI{m});
    end
    
end
    


% k-means clustering of the matrices of phase synch. (number of clusters
% using DBI) upperbound for DBI was set to nS 
[popidx{1},popCorr{1},~,~,~,DBI(1)] = mykmeans(cat(1,CCORSW{:}),nS,nS);
[popidx{2},popCorr{2},~,~,~,DBI(2)] = mykmeans(cat(1,COSDELPHI1{:}),nS,nS);

%%%%% for each simulation, match them with the population state --- or the
%%%%% ground truth
for i = 1:N
% k-means clustering of the matrices of phase synch.
[idx{1,i},Corr{1,i}] = mykmeans(CCORSW{i},nS,nS);
[idx{2,i},Corr{2,i}] = mykmeans(COSDELPHI1{i},nS,nS);
[~,~,mCorr{1,i},midx{1,i}] = matchstatesClutster(grdstate{1},Corr{1,i},idx{1,i},1);
[~,~,mCorr{2,i},midx{2,i}] = matchstatesClutster(grdstate{1},Corr{2,i},idx{2,i},1);
end

for i = 1:N
    for j = 1:nS 
        ccstate{j}(:,:,i) = mCorr{1,i}(:,:,j);
        crpstate{j}(:,:,i) = mCorr{2,i}(:,:,j);
    end
end

for i = 1:nS
    meanstate{i,1} = mean(ccstate{i},3);
    meanstate{i,2} = mean(crpstate{i},3);
end

count = 0;
placeholder = [1 2 3;7 8 9];
figure;
h = tight_subplot(3, nS, [.001 .001],[.01 .001],[.05 .05]);
labels = {'CIRC','CRP'};
set(h(4:6),'visible','off')
for j = 1:nS
    for i = 1:2
        count = count + 1;
        axes(h(placeholder(count)));gsplot(meanstate{j,i});
        axis square;
        set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w')
        if i == 1
         title(strcat(['State ' num2str(j)]));
        end
        if j == 1
            ylabel(labels{i},'interpreter','latex');
        end
    end
    
end
colorbar(h(8),'Position',[0.93 0.168 0.022 0.7]); caxis([-1 1]);


for i = 1:N
    MIDX1(:,i) = midx{1,i};
    MIDX2(:,i) = midx{2,i};
end

for i = 1:size(MIDX1,1)
    flag1(i) = sum( MIDX1(i,:) == grdidx(i))/N;
    flag2(i) = sum( MIDX2(i,:) == grdidx(i))/N;
end

for j = 1:nS
for i = 1:size(MIDX1,1)
    flags1(i,j) = sum( MIDX1(i,:) == j)/N;
    flags2(i,j) = sum( MIDX2(i,:) == j)/N;
end
end




figure;subplot(3,1,1);plot(grdidx,'r','LineWidth',2);ylim([0 4]);xlabel('t');
ylabel('State #');title('Ground truth state')
subplot(3,1,2);plot(flag1,'g','LineWidth',2);
title('Accuracy of correctly classifying the state');ylabel('Accuracy');xlabel('t')
subplot(3,1,3);bar(flags1,'stacked')
title('Stacked bar of the classification accuracy across time');
legend('State 1','State 2','State 3');
ylabel('State Accuracy Proportion');xlabel('t')




figure;imagesc(MIDX1');colormap winter;
imAlpha=ones(size(MIDX1'));
imAlpha(isnan(MIDX1'))=0;
imagesc(MIDX1','AlphaData',imAlpha);
set(gca,'color',1*[1 1 1]);
title('State transitions - na-MEMD-based-PS using tWPS: CIRC');
xlabel('time [s]');ylabel('realizations')
cb = colorbar();
% Set color labels (one for each row in RGB)
label = 1:nS; 
caxis([1,numel(label)])
cb.YTick = 1 : nS;
labelChar = label;
cb.TickLabels = labelChar(1:end-1);
cb.FontSize = 12; 

figure;imagesc(MIDX2');colormap winter
title('State transitions - na-MEMD-based-PS using IPS: CRP');
xlabel('time [s]');ylabel('realizations');
cb = colorbar();
% Set color labels (one for each row in RGB)
label = 1:nS; 
caxis([1,numel(label)])
cb.YTick = 1 : nS;
labelChar = label;
cb.TickLabels = labelChar(1:end-1);
cb.FontSize = 12; 


%% MVMD
indx = nchoosek(1:size(Data{1},2),2);
N=1000;
    
k = 4; % apriori knowledge from naMEMD
for m =1:N
    noise = mvnrnd([0 0 0],[1 0 0;0 1 0;0 0 1],length(t))';
    ex1 = noise(1,:);
    ex2 = noise(2,:);
    ex3 = noise(3,:);
    XN1 = x1 + ex1;
    XN2 = x2 + ex2;
    XN3 = x3 + ex3;
    Data{m} = [XN1;XN2;XN3]';
    % Using multivariate variational mode decomposition
    [u, u_hat, omega] = MVMD(Data{m}, 1000, 0.01, k, 1, 0, 1e-9);
    for i = 1:k
        mfreq{m}(1,i) = meanfreq(u(i,:,1),fs);
    end

    
     % Phase Synchronization analysis:
     [~,ind] = min(abs(mfreq{m}'-f));
    for i = 1:1:size(indx,1)
        dat{m} = [u(ind,:,indx(i,1));u(ind,:,indx(i,2))];
        H = hilbert(dat{m}');
        sigphase = angle(H);
        DELPHI{m} = sigphase(:,1)-sigphase(:,2);
        CCORSW{m}(:,i) = circularslidingwindow(sigphase(:,1),sigphase(:,2),winLen,'vonmises');%,'option','window','winsize',w(1));
        COSDELPHI1{m}(:,i) = cos(DELPHI{m});
    end

end
    



% k-means clustering of the matrices of phase synch. (number of clusters
% using DBI) upperbound for DBI was set to nS 
[popidx{1},popCorr{1},~,~,~,DBI(1)] = mykmeans(cat(1,CCORSW{:}),nS,nS);
[popidx{2},popCorr{2},~,~,~,DBI(2)] = mykmeans(cat(1,COSDELPHI1{:}),nS,nS);

%%%%% for each simulation, match them with the population state --- or the
%%%%% ground truth
for i = 1:N
% k-means clustering of the matrices of phase synch.
[idx{1,i},Corr{1,i}] = mykmeans(CCORSW{i},nS,nS);
[idx{2,i},Corr{2,i}] = mykmeans(COSDELPHI1{i},nS,nS);
[~,~,mCorr{1,i},midx{1,i}] = matchstatesClutster(grdstate{1},Corr{1,i},idx{1,i},1);
[~,~,mCorr{2,i},midx{2,i}] = matchstatesClutster(grdstate{1},Corr{2,i},idx{2,i},1);
end

for i = 1:N
    for j = 1:nS 
        ccstate{j}(:,:,i) = mCorr{1,i}(:,:,j);
        crpstate{j}(:,:,i) = mCorr{2,i}(:,:,j);
    end
end

for i = 1:nS
    meanstate{i,1} = mean(ccstate{i},3);
    meanstate{i,2} = mean(crpstate{i},3);
end

count = 0;
placeholder = [1 2 3;7 8 9];
figure;
h = tight_subplot(3, nS, [.001 .001],[.01 .001],[.05 .05]);
labels = {'CIRC','CRP'};
set(h(4:6),'visible','off')
for j = 1:nS
    for i = 1:2
        count = count + 1;
        axes(h(placeholder(count)));gsplot(meanstate{j,i});
        axis square;
        set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w')
        if i == 1
         title(strcat(['State ' num2str(j)]));
        end
        if j == 1
            ylabel(labels{i},'interpreter','latex');
        end
    end
    
end
colorbar(h(8),'Position',[0.93 0.168 0.022 0.7]); caxis([-1 1]);


for i = 1:N
    MIDX1(:,i) = midx{1,i};
    MIDX2(:,i) = midx{2,i};
end

for i = 1:size(MIDX1,1)
    flag1(i) = sum( MIDX1(i,:) == grdidx(i))/N;
    flag2(i) = sum( MIDX2(i,:) == grdidx(i))/N;
end

for j = 1:nS
for i = 1:size(MIDX1,1)
    flags1(i,j) = sum( MIDX1(i,:) == j)/N;
    flags2(i,j) = sum( MIDX2(i,:) == j)/N;
end
end




figure;subplot(3,1,1);plot(grdidx,'r','LineWidth',2);ylim([0 4]);xlabel('t');
ylabel('State #');title('Ground truth state')
subplot(3,1,2);plot(flag1,'g','LineWidth',2);
title('Accuracy of correctly classifying the state');ylabel('Accuracy');xlabel('t')
subplot(3,1,3);bar(flags1,'stacked')
title('Stacked bar of the classification accuracy across time');
legend('State 1','State 2','State 3');
ylabel('State Accuracy Proportion');xlabel('t')




figure;imagesc(MIDX1');colormap winter;
imAlpha=ones(size(MIDX1'));
imAlpha(isnan(MIDX1'))=0;
imagesc(MIDX1','AlphaData',imAlpha);
set(gca,'color',1*[1 1 1]);
title('State transitions - MVMD-based-PS using tWPS: CIRC');
xlabel('time [s]');ylabel('realizations')
cb = colorbar();
% Set color labels (one for each row in RGB)
label = 1:nS; 
caxis([1,numel(label)])
cb.YTick = 1 : nS;
labelChar = label;
cb.TickLabels = labelChar(1:nS);
cb.FontSize = 12; 

figure;imagesc(MIDX2');colormap winter
title('State transitions - MVMD-based-PS using IPS: CRP');
xlabel('time [s]');ylabel('realizations');
cb = colorbar();
% Set color labels (one for each row in RGB)
label = 1:nS; 
caxis([1,numel(label)])
cb.YTick = 1 : nS;
labelChar = label;
cb.TickLabels = labelChar(1:nS);
cb.FontSize = 12; 