% Code for Figure 2
% Bias of the estimator for a Noiseless Signal

N = 8; %8 - Point DFT
n = 0:1:N-1;
xn = 5*(exp(1i*(pi*0.65)*n)); %Same signal used in Figure 1
DFT = fft(xn,N);
[Max,kp] = max(abs(DFT)); %Getting index of peak frequency of the signal

%Initializing arrays to store Bias Values
biasJacobsen = zeros(1,19);
biasCandan = zeros(1,19);
biasParabolic = zeros(1,19);
biasMacleod = zeros(1,19);
biasQuinn = zeros(1,19);

%ind is used to iterate through the Bias Arrays
ind = 1;

%{
Method Used:

We will vary delta for a fixed range, and for each delta we will generate a signal
and compute its DFT(usign fft()). 
Now using the DFT samples, we can compute R[kpeak-1], R[kpeak] and R[kpeak+1]
Using these three DFT samples, we can estimate delta usign the different
estimators(Macleod, Quinn, Jacobsen, Candan(proposed))

We can store absolute difference in estimate delta and actual delta in the
bias array
%}


for d=0.025:0.025:0.475 % Varying delta 
    w = 2*pi*(kp+d)/N;
    r = 5*(exp(1i*w*n));
    R = fft(r,N);
    
    %We add one to the index obtained from max(), because MATLAB has
    %1-based indexing but DFT has 0-based indexing
    Rkp = R(kp+1); %R(kpeak)
    RkpM = R(kp); %R(kpeak-1)
    RkpP = R(kp+2); %R(kpeak+1)
    
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
    
    %Calculting Bias = E[y^] - y
    biasCandan(ind) = abs(deltaCandan-d); %  
    biasJacobsen(ind) = abs(deltaJacobsen-d);  
    biasParabolic(ind) = abs(deltaParabolic-d);  
    biasMacleod(ind) = abs(deltaMacleod-d);
    biasQuinn(ind) = abs(dQuinn-d);
    
    %incrementing index of bias array
    ind = ind +1;
end

d=0.025:0.025:0.475;  
hold on;
plot(d, biasJacobsen, '-d','color','#666699','LineWidth',1.5);
plot(d, biasParabolic,'color','#339933','LineWidth',1.5);
plot(d, biasMacleod , '--s','color','#cc0000','LineWidth',1.5);
plot(d, biasQuinn , '--x','color','#0033cc','LineWidth',1.5);
plot(d, biasCandan, '-o','color','#993399','LineWidth',1.5);

hold off;
set(gca, 'YScale', 'log')
ylim([0.0000001 1]);
xlim([0 0.5]);
legend('Jacobsen','Parabolic','MacLeod','Quinn','Proposed');
ylabel('Bias'); xlabel('\delta');
