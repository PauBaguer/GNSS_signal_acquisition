%EXERCISE 2: GPS L1 CA signal analysis
% Luis Esteve NACC 2023
% luis.esteve@upc.edu
% Javier Arribas NACC 2018
% jarribas@cttc.es

clear all;
close all;
addpath('../../libs');
%% DEFINE THE SIGNAL PARAMETERS
%Define the sampling frequency [Hz]
Fs=10e6;
%Define the digitized signal duration [s]
Tsig=0.001; 
%Define the carrier frequency
Fc=Fs/4;
%define chip frequency
Rc_L1_CA=1.023e6;

% **** QUESTION 1 ****: 
% Compute the number of chips that fits in Tsig and store it in num_chips
%num_chips=?

%generate the time vector
Ts=1/Fs;
t=0:Ts:(Tsig-Ts);
Nsamples=length(t);

%% generate the transmitted satellite baseband signal
NSATS=4;
 for numsat=1:1:NSATS
     % Generate GPS L1 CA signal using our custom generator function.
     % **** QUESTION 2 ****:
     % Open the source code of digitGPS_L1_CA.m and codegen.m. Explain what
     % is the purpose of each function.
    s_sat_tx(numsat,:)=digitGPS_L1_CA(num_chips,Fs,0,numsat,1); % PRN C_E1_B muestreada a Rc_E1_B (directamente chip rate)
 end
 % store the baseband signal of the first satellite
 s_bb_tx=s_sat_tx(1,:);
 
 % **** QUESTION 3 ****: 
% Perform the FFT of s_bb_tx. 
% Check MATLAB FFT help by typing help FFT
% S_bb_tx=?

%compute the frequency axis in Hz
f = Fs/2*linspace(-1,1,Nsamples);
%reorder the FFT results and compute the signal power spectrum in dBW/Hz
NFFT=Nsamples;
S_bb_tx_spectrum_db=10*log10((abs([S_bb_tx(floor((NFFT/2))+1:1:end) S_bb_tx(1:1:(NFFT/2))]).^2)/Nsamples);

%plot the spectrum
% **** QUESTION 4 ****: 
% Capture a figure showing the spectrum of the GPS L1 CA signal. 
% Adjust the zoom to be able to Identify the main lobes and the secondary lobes of the modulation. 
% Measure the distance in Hz from the center frequency to the first spectrum
% null of the GPS L1 CA modulation, each side.
figure;
plot(f,S_bb_tx_spectrum_db);
title('Transmitted baseband signal');
xlabel('Hz');
ylabel('dBW/Hz');
 
%% Generate a RF signal carying the baseband signal
% notice the relation of the signal power of a sinusoid:
% power(s(t))= [s(t)]^2 -> a*cos(w)^2=(a^2/2)+(a^2/2)*cos(2w)
% average_power=avr((a^2/2)+(a^2/2)*cos(2w))=(a^2/2)

%Define the digitized signal amplitude
Psig=1;


% **** QUESTION 5 ****: 
% Write here the equation that produces an RF signal, modulating a
% carrier frequency Fc with the baseband signal s_bb_tx. Store it in 
% a vector named s_rf_tx. The obtained signal must have a signal power
% equal to Psig.

%Define the digitized signal amplitude
%Asig=?

%Compute the RF transmited signal
%s_rf_tx=?

% Perform the FFT
S_rf_tx=fft(s_rf_tx);
%compute the frequency axis in Hz
f = Fs/2*linspace(-1,1,Nsamples);
%reorder the FFT results and compute the signal power spectrum in dBW/Hz
NFFT=Nsamples;
S_rf_tx_spectrum_db=10*log10((abs([S_rf_tx(floor((NFFT/2))+1:1:end) S_rf_tx(1:1:(NFFT/2))]).^2)/Nsamples);
%plot the spectrum
figure;
plot(f,S_rf_tx_spectrum_db);
title('Transmitted RF signal');
xlabel('Hz');
ylabel('dBW/Hz');

% **** QUESTION 6 ****: 
% Capture the figure at maximum zoom out.

%% simulate the received signal and add noise

% Doppler frequency
Fd=2e3;

% **** QUESTION 7 ****: 
% Write here the equation that produces an RF received signal, modulating a
% carrier frequency Fc with the baseband signal s_bb_tx, and affected by a
% Doppler shift Fd. Store it in a vector named s_rf_rx. The obtained signal
% must have a signal power equal to Psig.
% s_rf_rx=?

SNR_dB=3;
SNR_lin=10^(SNR_dB/10);
An=sqrt(Psig/SNR_lin);
noise=An*randn(1,Nsamples);
BB_BW_Hz=Fs/2.01;
[b,a] = butter(5,BB_BW_Hz/(Fs/2),'low');
noise=filtfilt(b,a,noise);
y_rf_rx=s_rf_tx+noise;

% Perform the FFT
Y_rf_tx=fft(y_rf_rx);
%compute the frequency axis in Hz
f = Fs/2*linspace(-1,1,Nsamples);
%reorder the FFT results and compute the signal power spectrum in dBW/Hz
NFFT=Nsamples;
Y_rf_tx_spectrum_db=10*log10((abs([Y_rf_tx(floor((NFFT/2))+1:1:end) Y_rf_tx(1:1:(NFFT/2))]).^2)/Nsamples);
%plot the spectrum
figure;
plot(f,Y_rf_tx_spectrum_db);
title('Received RF signal');
xlabel('Hz');
ylabel('dBW/Hz');

% **** QUESTION 8 ****: 
% Capture the figure at maximum zoom out. Explain the effects of noise.

% **** QUESTION 9 ****: 
% Perform a zoom in and measure the approximated center frequency 
% of the received signal. It is possible?


%% SNR estimation
% **** QUESTION 10 ****: 
% Compute the SNR and the CNO in dB
% HINT: Consider RF_BW_Hz=Fs/2.
% HINT: Remember that the signal and noise power can be estimated by using the
% autocorrelation estimation Rxx=(1/K)*(x*x'), where x is a vector
% containing the signal or noise samples.

% SNR_estim_dB=?

disp(['Estimated SNR is ' num2str(SNR_estim_dB) ' vs. desired SNR: ' num2str(SNR_dB) ' [dB]']);
% CN0_estim_dB=?

disp(['Estimated CN0 is ' num2str(CN0_estim_dB) ' [dB]']);

%% Downconvert the RF signal to baseband
% design the filters
% Low pass filter BW.
% **** QUESTION 11 ****: 
% Compute the maximum allowable baseband BW that prevents the signal
% aliasing after the donwconversion to baseband, considering an RF bandwith of Fs/2.
% BB_BW_Hz=?

[b,a] = butter(3,BB_BW_Hz/(Fs/2),'low');

% **** QUESTION 12 ****: 
% Perform the downconversion and IQ demodulation, considering
% the expected carrier frequency Fc. Store the baseband signal in s_bb_rx
% HINT: The Low Pass Filter (LPF) can be applied by using
% y=filtfilt(b,a,x), where x is the input signal vector and y is the filtered
% output signal vector.

% s_bb_rx_I_filt=?
% s_bb_rx_Q_filt=?

%ensamble the recovered baseband signal (complex envelope)
s_bb_rx=s_bb_rx_I_filt+1j*s_bb_rx_Q_filt;

% Perform the FFT
S_bb_rx=fft(s_bb_rx);
%plot the spectrum magnitude
%compute the frequency axis in Hz
f = Fs/2*linspace(-1,1,Nsamples);
%reorder the FFT results and compute the signal power spectrum in dB
NFFT=Nsamples;
S_bb_rx_spectrum_db=10*log10((abs([S_bb_rx(floor((NFFT/2))+1:1:end) S_bb_rx(1:1:(NFFT/2))]).^2)/Nsamples);
figure;

% **** QUESTION 13 ****: 
% Capture the received baseband signal spectrum figure. Can you see the
% main lobe?
plot(f,S_bb_rx_spectrum_db);
title('Received baseband signal');
xlabel('Hz');
ylabel('dBW/Hz');

%% Real data
% **** QUESTION 14 ****:
% Compute the SNR in dB of a satellite signal with a nominal CN0 on the Earth
% surface of 45 dB-Hz. Set the obtained value in SNR_dB in the noise
% generation section of this code and regenerate the plots. What can you
% see in the new plots?



