clc;
close all;
clear all;
rng('shuffle');
% Reinitialize the random number generator used by RAND ,randi, and RANDN
SNRdB = [1,2,3];
m = 5;
numBlocks = 10000;%no.of times sensing
ssgam = logspace(0.005,1.5,50);%generates 50 points between 10^0.005 and 10^1 in log scale,threshold samples
%empty vectors with same size as threshold vectors vector.
ssDet = zeros(length(SNRdB),length(ssgam));
%this code is for unit variance cond
ssFA = zeros(length(SNRdB),length(ssgam));
SNR = zeros(length(SNRdB),length(ssgam));
for i=1:length(SNRdB)
for L = 1:numBlocks
	pisym=2*randi([0,1],[m,1])-1;%m no.of binary phase shifted symbols(X).
	h=1/sqrt(2)*(randn+j*randn);%normalized fading coeff(H).
	ChNoise = 1/sqrt(2)*(randn(m,1)+j*randn(m,1));%white gaussian noise(N)
	for K = 1 : length(ssgam)
	SNR(i,K) = 10^(SNRdB(i)/10);%converting SNR to normal form from DB.
	TxVec = sqrt(SNR(i,K)*pisym);%transmitted vector from CR
	RxVec_H1 = h*TxVec + ChNoise;%rxed vector at energy detector for H1 case
    ssStat_H1 = abs(TxVec'/norm(TxVec)*RxVec_H1)^2;%calculating energy of rxed signal for H1 case
    ssDet(i,K) = ssDet(i,K) + (ssStat_H1>=ssgam(K));%total no.of detections for different Thresholds
    RxVec_H0 =  ChNoise;%rxed vector at energy detector for H0 case
    ssStat_H0 = abs(TxVec'/norm(TxVec)*RxVec_H0)^2;%calculating energy of rxed signal for H0 case
	ssFA(i,K) = ssFA(i,K) + (ssStat_H0>=ssgam(K));%total no.of false alarms  for different Thresholds
	end
end
ssDet(i,:) = ssDet(i,:)/(numBlocks);%probability of detection for particular SNR
ssFA(i,:) = ssFA(i,:)/(numBlocks);%probability of false alarm
end
 plot(ssFA(1,:) , ssDet(1,:), 'b-', 'linewidth',2.0);
 hold on;
 plot(exp(-ssgam),exp(-ssgam./(m*SNR(1,:)+1)),'b o', 'linewidth',2.0);
 hold on;
plot(ssFA(2,:) , ssDet(2,:),'r-', 'linewidth',2.0);
 hold on;
 plot(exp(-ssgam),exp(-ssgam./(m*SNR(2,:)+1)),'r o', 'linewidth',2.0);
 hold on;
 plot(ssFA(3,:) , ssDet(3,:), 'g-', 'linewidth',2.0);
 hold on;
 plot(exp(-ssgam),exp(-ssgam./(m*SNR(3,:)+1)),'g o', 'linewidth',2.0);
axis tight;
grid on;
legend('sim with SNR=1','Theory with SNR=1','sim with SNR=2','Theory with SNR=2','sim with SNR=3','Theory with SNR=3');
xlabel('P_{FA}');
ylabel('P_D');
title('Receiver operating charactrstic');
