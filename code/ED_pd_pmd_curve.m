clear all; 
close all;
clc;
n=10; %Number of data symbols
Tsym=8; %Symbol time interms of sample time or oversampling rate equivalently
SNRstep=5;
SNR_dB=-20:SNRstep:20;%takes SNRdB values from [-20 -15 -10 -5 0 5 10 15 20]
rng('default');%set the random generator seed to default.
BinData=round(rand(1,n));%generates n no.of binary numbers randomly
data=2*BinData-1; %BPSK data Unipolar 2 Bipolar
bpsk=reshape(repmat(data,1,Tsym)',n*Tsym,1); %BPSK signal with required duty cycle
L=length(bpsk);
PDc2=0;
PDc1=0;
PFc=0;
PFc1=0;
PFc2=0;
ssDet=0;
ssFa=0;
for i=1:length(SNR_dB)
    SNR = 10^(SNR_dB(i)/10); %SNR dB to linear scale
    h=1/sqrt(2)*(randn+j*randn);
    Esym=sum(abs(bpsk).^2)/(L); %Calculate actual symbol energy
    N0=Esym/SNR; %Find the noise spectral density
    noiseSigma = sqrt(N0);%Standard deviation for AWGN Noise when x is real
    noise(:,i) = noiseSigma*randn(1,L);%computed noise

    receivedx(:,i)=bpsk+noise(:,i) ;%received signal with only AWGN or ip to matcehd filter
    receivedfliped=flipud(receivedx(:,i));%flipped received signal or matched filter responsse

    impRes = [0.5 ones(1,6) 0.5]; %Averaging Filter -&gt; u[n]-u[n-Tsamp]
    yy1(:,i)=conv(receivedx(:,i),impRes,'full');%matched filter output

    xp(:,i)=receivedfliped;
    thresh(:,i)=noise(:,i).*xp(:,i);
    yy2(:,i)=filter(xp(:,i),1,receivedx(:,i));
    TxVec(:,i) = sqrt(SNR*bpsk);%transmitted vector from CR
	RxVec_H1(:,i) = h*TxVec(:,i) + noise(:,i);%rxed vector at energy detector for H1 case
    ssStat_H1 (:,i)= abs((TxVec(:,i)'/norm(TxVec(:,i))).*RxVec_H1(:,i)').^2;%calculating energy of rxed signal for H1 case
    RxVec_H0 (:,i)= noise(:,i);%rxed vector at energy detector for H0 case
    ssStat_H0 (:,i)= abs((TxVec(:,i)'/norm(TxVec(:,i))).*RxVec_H0(:,i)').^2;%calculating energy of rxed signal for H0 case
	

    for j = 1:L
        if(ssStat_H1(j,i)>thresh(j,i))
            ssDet = ssDet+1;
        end
        if(ssStat_H0(j,i)>thresh(j,i))
            ssFa= ssFa +1;
        end
        if(yy1(j,i)>thresh(j,i))
            PDc1=PDc1+1;
        end
        if(yy2(j,i)>thresh(j,i))
            PDc2=PDc2+1;
        end
        if(noise(j,i)>thresh(j,i))
            PFc=PFc+1;
        end
         if(yy1(j,i)<thresh(j,i))
            PFc1=PFc1+1;
        end
        if(yy2(j,i)<thresh(j,i))
            PFc2=PFc2+1;
        end
    end
    
     PD2(:,i)=PDc2/(L*n);
     PFinv(:,i)=(PFc/(L*n));
     PF(:,i)=1-PFinv(:,i);
     Det(:,i)=ssDet/(L*n);
      Fa(:,i)=ssFa/(L*n);
   
end
pmd1=1-PD2;
pmd2=1-Det;
plot(SNR_dB,Det,'-bo',...
    'LineWidth',2,'MarkerEdgeColor','y',...
    'MarkerFaceColor','g','MarkerSize',5);
hold on
plot(SNR_dB,pmd2,'-ro',...
    'LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',5);
legend('detection curve','missed detection curve ')
ylabel('Probability of detection and missed detection ');
xlabel('SNR_{dB}');
title('Probability of detection and missed detection  curves in energy detector');
