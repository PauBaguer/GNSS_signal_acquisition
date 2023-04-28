%EXERCISE 1: front-end simulator and signal analysis
% Luis Esteve NACC 2023
% luis.esteve@upc.edu
% Javier Arribas NACC 2018
% jarribas@cttc.es

clear all;
close all;
addpath('../../libs');
%% DEFINE THE SIGNAL PARAMETERS
%Define the sampling frequency [Hz]
% Fs>2Fc
Fs=2e6;
%Define the digitized signal duration [s]
Tsig=0.1; 
%Define the carrier frequency
%In order to not overload the MATLAB with high sampling frequency, 
% we use here a low carrier frequency.
Fc=500e3;

%generate the time vector
Ts=1/Fs;
t=0:Ts:(Tsig-Ts);
Nsamples=length(t);

%% Generate a BPSK baseband signal
initial_phase_rad=pi/4;
%BPSK symbol rate
f_symbol_hz=1e3;
%produce random symbols
bpsk_symbols=2*(rand(1,floor(f_symbol_hz*Tsig))>0.5)-1;
bpsk_symbols_resampled=resample(bpsk_symbols,Fs,f_symbol_hz,0);
s_bb_tx=bpsk_symbols_resampled*exp(1j*initial_phase_rad);

% Perform the FFT
% **** QUESTION 1 ****: 
% Perform the FFT of s_bb_tx.  
% Check MATLAB FFT help by typing help FFT
%S_bb_tx=?


%compute the frequency axis in Hz
f = Fs/2*linspace(-1,1,Nsamples);
%reorder the FFT results and compute the signal power spectrum in dBW/Hz
NFFT=Nsamples;
S_bb_tx_spectrum_db=10*log10((abs([S_bb_tx(floor((NFFT/2))+1:1:end) S_bb_tx(1:1:(NFFT/2))]).^2)/Nsamples);

%plot the spectrum
% **** QUESTION 2 ****: 
% Capture a figure showing the spectrum of the BPSK signal from -3kHz to 3kHz. 
% Identify the main lobes and the secondary lobes of the modulation. 
% Which parameter of the BPSK modulation has an effect on the lobe spacing?
% Check the effect by variying such parameter and perform the plot again.
figure;
plot(f,S_bb_tx_spectrum_db);
title('Transmitted baseband signal');
xlabel('Hz');
ylabel('dBW/Hz');


%% Generate a RF signal carrying the baseband signal
% notice the relation of the signal power of a sinusoid:
% power(s(t))= [s(t)]^2 -> a*cos(w)^2=(a^2/2)+(a^2/2)*cos(2w)
% average_power=avr((a^2/2)+(a^2/2)*cos(2w))=(a^2/2)

%Define the RF signal power
Psig=4;


% **** QUESTION 3 ****: 
% Write here the equation that produces an RF signal, modulating a
% carrier frequency Fc with the baseband signal s_bb_tx. Store it in 
% a vector named s_rf_tx. The obtained signal must have a signal power
% equal to Psig.

%Compute the signal amplitude
%Asig=?

%Compute the RF signal
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

% **** QUESTION 4 ****: 
% Capture the figure at maximum zoom out. Explain why it is symetrical.

%% simulate the received signal and add noise

%Define a Doppler shift
Fd=2e3;

% **** QUESTION 5 ****: 
% Write here the equation that produces an RF received signal, modulating a
% carrier frequency Fc with the baseband signal s_bb_tx, and affected by a
% Doppler shift Fd. Store it in a vector named s_rf_rx. The obtained signal
% must have a signal power equal to Psig.

%s_rf_rx=?

%add noise
SNR_dB=3;
SNR_lin=10^(SNR_dB/10);

% **** QUESTION 6 ****: 
% Compute the noise amplitude (linear)
%An=?


noise=An*randn(1,Nsamples);
BB_BW_Hz=Fs/2.1;
[b,a] = butter(5,BB_BW_Hz/(Fs/2),'low');
noise=filtfilt(b,a,noise);
y_rf_rx=s_rf_rx+noise;

% Perform the FFT
Y_rf_rx=fft(y_rf_rx);
%compute the frequency axis in Hz
f = Fs/2*linspace(-1,1,Nsamples);
%reorder the FFT results and compute the signal power spectrum in dBW/Hz
NFFT=Nsamples;
Y_rf_rx_spectrum_db=10*log10((abs([Y_rf_rx(floor((NFFT/2))+1:1:end) Y_rf_rx(1:1:(NFFT/2))]).^2)/Nsamples);
%plot the spectrum
figure;
plot(f,Y_rf_rx_spectrum_db);
title('Received RF signal (Corrección Luis).');
xlabel('Hz');
ylabel('dB/Hz');

% **** QUESTION 7 ****: 
% Capture the figure at maximum zoom out. Explain the effects of noise.

% **** QUESTION 8 ****: 
% Perform a zoom in and measure the approximated center frequency 
% of the received signal.

%% SNR estimation
% **** QUESTION 9 ****: 
% Compute the SNR and the CNO in dB
% HINT: Consider RF_BW_Hz=Fs/2.
% HINT: Remember that the signal and noise power can be estimated by using the
% autocorrelation estimation Rxx=(1/K)*(x*x'), where x is a vector
% containing the signal or noise samples.
%SNR_estim_dB=?

disp(['Estimated SNR is ' num2str(SNR_estim_dB) ' vs. desired SNR: ' num2str(SNR_dB) ' [dB]']);

%CN0_estim_dB=?

disp(['Estimated CN0 is ' num2str(CN0_estim_dB) ' [dB]']);

%% Downconvert the RF signal to baseband
% design the LPF filters
BB_BW_Hz=Fc/2;
[b,a] = butter(5,BB_BW_Hz/(Fs/2),'low');

% **** QUESTION 10 ****: 
% Perform the downconversion and IQ demodulation, considering
% the expected carrier frequency Fc. Store the baseband signal in s_bb_rx
% HINT: The Low Pass Filter (LPF) can be applied by using
% y=filtfilt(b,a,x), where x is the input signal vector and y is the filtered
% output signal vector.

%s_bb_rx_I_filt=?
%s_bb_rx_Q_filt=?

%ensamble the recovered baseband signal (complex envelope)
s_bb_rx=s_bb_rx_I_filt+1j*s_bb_rx_Q_filt;

% Perform the FFT
S_bb_rx=fft(s_bb_rx);
%plot the spectrum magnitude
%compute the frequency axis in Hz
f = Fs/2*linspace(-1,1,Nsamples);
%reorder the FFT results and compute the signal power spectrum in dBW/Hz
NFFT=Nsamples;
S_bb_rx_spectrum_db=10*log10((abs([S_bb_rx(floor((NFFT/2))+1:1:end) S_bb_rx(1:1:(NFFT/2))]).^2)/Nsamples);
figure;

% **** QUESTION 11 ****: 
% Capture the received baseband signal spectrum figure. Estimate
% the approximated LPF filter bandwidth and estimate the approximated BPSK signal bandwidth

plot(f,S_bb_rx_spectrum_db);
title('Received baseband signal (Corrección Luis)');
xlabel('Hz');
ylabel('dBW/Hz');




