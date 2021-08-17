%Code for Figure 1 - To demonstrate difference in frequency of signal from
%DFT vs DTFT

clear; clc;

N = 16; %16 Point DFT
n = 0:1:N-1;
r = 5*(exp(1i*(pi*0.65)*n));

w = 0:0.01:2*pi;
dtft = DT_Fourier(r,0,w);
R = fft(r,N);
hold on;
stem(n,abs(R),'LineWidth',1.5);
plot(w*(N/(2*pi)),abs(dtft),'LineWidth',1.5);

[Max,kp] = max(abs(R));

x = [kp-2 kp-1 kp];
y = [abs(R(kp-1)) abs(R(kp)) abs(R(kp+1))];

%Labelling
labels = {'|R[k_{p}-1]|','|R[k_{p}]|','|R[k_{p}+1]|'};
labels2 = {'k_{p}-1','k_{p}','k_{p}+1'};
text(x,y,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
text(x,[-2 -2 -2],labels2,'VerticalAlignment','top','HorizontalAlignment','center')

[Max,peak] = max(dtft);
s = size(w);
sizew = s(2);

stem(peak*N/sizew,abs(Max),'LineWidth',1.5);
text(peak*N/sizew,-0.5,'k_{p}+\delta','VerticalAlignment','top','HorizontalAlignment','center')
xlabel('w'); ylabel('|X(e^{jw})|');
legend('DFT', 'DTFT');

