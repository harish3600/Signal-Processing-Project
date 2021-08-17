clear; clc;

N = 8; %N Point DFT
n = 0:1:N-1; 
xn = 1*(exp(1i*(pi*0.65)*n)); %Signal x[n]
signalPower = sum(abs(xn.^2)); %Signal Power of x[n]

SNR = zeros(1,80); %Vector to store SNR values

%Initialsing vectors to store RMSE values for different Estimators
RMSEJacobsen = zeros(1,80);
RMSECandan = zeros(1,80);
RMSEParabolic = zeros(1,80);
RMSEQuinn = zeros(1,80);
RMSEMacleod = zeros(1,80);


d = 0.25; %Given delta = 0.25
I = 1; %Variable to iterate the RMSE and SNR array
T = 1000;

%Calculating noise power for varying SNR
for SN=0:3:81 %Varying SNR
  Ps = signalPower; %Fixing the signal power
  biasjsum = 0;
  biascsum = 0;
  biaspsum = 0;
  biasqsum = 0;
  biasmsum = 0;
  
  for int=0:T %Running trials by varying phase for signals 
    Pn = Ps*(10^(-SN/10)); %Calculating power of noise from SNR and Power of Signal
    nzv=sqrt(Pn)*randn(1,N); %Generating a noise signal 
    
    %Get kpeak of the signal = xn + noise
    sig = xn+nzv; %Adding noise to the input signal
    DFT = fft(sig,N);
    [Max,kp] = max(abs(DFT)); %DFT of Noisy Signal
    
    w = 2*pi*(kp+d)/N; %Generating a new signal with frequency such that delta = 0.25     
    r = 1*(exp(1i*w*n));
    rn = r + nzv; %r[n] = Ae + w[n]
  
    R = fft(rn,N);
    
    %Using mod, as DFT if N perioic
    %This avoids out-of-bound error when kp is the last index in R
    
    Rkp = R(mod(kp,N)+1);
    RkpM = R(kp);
    RkpP = R(mod(kp+1,N)+1);
  
    %We add one to the index obtained from max(), because MATLAB has
    %1-based indexing but DFT has 0-based indexing
    
    %Estimating delta by using different estimators formula
    
    deltaJacobsen = real((RkpM-RkpP)/(2*Rkp-RkpM-RkpP));
    
    biasCorrectionfactor = (tan(pi/N)*N/pi);
    deltaCandan = biasCorrectionfactor*real((RkpM-RkpP)/(2*Rkp-RkpM-RkpP));
    
    deltaParabolic = ((abs(RkpP)-abs(RkpM))/(4*abs(Rkp)- 2*abs(RkpM) - 2*abs(RkpP)));
    
    dMac = real(RkpM*conj(Rkp)-RkpP*conj(Rkp))/real(2*abs(Rkp).^2 + RkpM*conj(Rkp) + RkpP*conj(Rkp));
    deltaMacleod = (sqrt(1+8*dMac^2)-1)/(4*dMac);
    
    %Quinn:
    a1 = real(RkpM/Rkp);
    a2 = real(RkpP/Rkp);
    dq1 = a1/(1-a1);
    dq2 = a2/(1-a2);
    if dq1>0 && dq2>0
        dQuinn = dq2;
    else
        dQuinn = dq1;
    end
    
    
    biasjsum = biasjsum + abs(deltaJacobsen-d)^2;
    biascsum = biascsum + abs(deltaCandan-d)^2;
    biaspsum = biaspsum + abs(deltaParabolic-d)^2;
    biasmsum = biasmsum + abs(deltaMacleod-d)^2;
    biasqsum = biasqsum + abs(dQuinn-d)^2;
  end
  RMSEJacobsen(I) = sqrt(biasjsum/T);
  RMSECandan(I) = sqrt(biascsum/T);
  RMSEMacleod(I) = sqrt(biasmsum/T);
  RMSEQuinn(I) = sqrt(biasqsum/T);
  RMSEParabolic(I) = sqrt(biaspsum/T);

  SNR(I) = SN;
  I = I+1;
  
end				

hold on;

plot(SNR, RMSEJacobsen, '-d','color','#666699','LineWidth',1.5);
plot(SNR, RMSEParabolic,'color','#339933','LineWidth',1.5);
plot(SNR, RMSEMacleod , '--s','color','#cc0000','LineWidth',1.5);
plot(SNR, RMSEQuinn , '--x','color','#0033cc','LineWidth',1.5);
plot(SNR, RMSECandan, '-o','color','#993399','LineWidth',1.5);

set(gca, 'YScale', 'log')
ylim([0.00001 100]);
xlim([0 80]);
xlabel('Input SNR');
ylabel('RMSE');
title('RMSE N=8 del=0.25');
legend('Jacobsen','Parabolic','MacLeod','Quinn','Proposed');
