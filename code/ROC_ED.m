clc;
close all;
clear all;
rng('shuffle');
% Reinitialize the random number generator used by RAND ,randi, and RANDN
SNRdB =2;
m = 5;
numBlocks = 10000;%no.of times sensing
ssgam = logspace(0.005,2,50);%generates 50 points between 10^0.005 and 10^1 in log scale,threshold samples
ssDet = zeros(size(ssgam));%empty vectors with same size as threshold vectors vector.
ssFA = zeros(size(ssgam));%this code is for unit variance cond
SNR = zeros(size(ssgam));
for L = 1:numBlocks
	pisym=2*randi([0,1],[m,1])-1;%m no.of binary phase shifted symbols(X).
	h=1/sqrt(2)*(randn+j*randn);%normalized fading coeff(H).
	ChNoise = 1/sqrt(2)*(randn(m,1)+j*randn(m,1));%white gaussian noise(N)
	for K = 1 : length(ssgam)
	SNR(K) = 10^(SNRdB/10);%converting SNR to normal form from DB.
	TxVec = sqrt(SNR(K)*pisym);%transmitted vector from CR
	RxVec_H1 = h*TxVec + ChNoise;%rxed vector at energy detector for H1 case
    ssStat_H1 = abs(TxVec'/norm(TxVec)*RxVec_H1)^2;%calculating energy of rxed signal for H1 case
    ssDet(K) = ssDet(K) + (ssStat_H1>=ssgam(K));%total no.of detections for different Thresholds
    RxVec_H0 =  ChNoise;%rxed vector at energy detector for H0 case
    ssStat_H0 = abs(TxVec'/norm(TxVec)*RxVec_H0)^2;%calculating energy of rxed signal for H0 case
	ssFA(K) = ssFA(K) + (ssStat_H0>=ssgam(K));%total no.of false alarms  for different Thresholds
	end
end
ssDet = ssDet/(numBlocks);%probability of detection
ssFA = ssFA/(numBlocks);%probability of false alarm
plot(ssFA , ssDet, 'b-', 'linewidth',2.0);
hold on;
plot(exp(-ssgam),exp(-ssgam./(m*SNR+1)),'r o', 'linewidth',2.0);
axis tight;
grid on;
xlabel('P_{FA}');
ylabel('P_D');
title('Receiver Operating Characteristic in energy detector');
