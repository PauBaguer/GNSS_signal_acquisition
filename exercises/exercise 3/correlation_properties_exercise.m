%EXERCISE 3: GPS L1 CA correlation properties
% Luis Esteve NACC 2021
% luis.esteve@upc.edu
% Javier Arribas NACC 2018
% jarribas@cttc.es

clear all;
close all;
addpath('../../libs');
%% DEFINE THE SIGNAL PARAMETERS
%Define the sampling frequency [Hz]
Fs=8e6;
%Define the digitized signal duration [s]
T_prn=1e-3;
Tsig=T_prn; 
%define chip frequency
Rc_L1_CA=1.023e6;
% **** QUESTION 1 ****:
% Compute the number of chips that fits in Tsig and store it in num_chips
% Considering that the GPS L1 CA PRN sequence has 1023 chips and a chip rate of 1.023 MChips/s, compute the
% duration of the GPS L1 CA complete sequence.
% num_chips=?
pause;

%generate the time vector
Ts=1/Fs;
t=0:Ts:(Tsig-Ts);
Nsamples=length(t);

%% generate the transmitted satellite baseband signal
NSATS=4;
 %satellite PRN 1
 for numsat=1:1:NSATS
    s_bb_sat(numsat,:)=digitGPS_L1_CA(num_chips,Fs,0,numsat,1); % PRN C_E1_B muestreada a Rc_E1_B (directamente chip rate)
 end
 
%% GPS L1 CA autocorrelation Rdd

d=s_bb_sat(1,:);
shift_samples=-floor(Nsamples/2):1:floor(Nsamples/2);

% **** QUESTION 2 ****:
% compute the circular autocorrelation estimation of the satellite signal stored in the vector
% d for all the possible delays, using the formula
% Rdd(tau)=(1/K)*(d*d(tau)'), where tau is the delay in samples. All the
% possible delay values are stored in the vector shift_samples.
% HINT: in order to perform a circular shift of a vector, use the built-in
% MATLAB function circshift. Use the MATLAB help to learn how to use
% circshift

for n=xxxxx
%     Rdd(n)=?;
end
pause;

% **** QUESTION 3 ****:
% Analize the following two figures. Check the X axis units. Measure the
% autocorrelation lag (or code phase shift) in CHIPS UNITS that produces a reduction of
% 50% of the autocorrelation gain. Convert the result to samples.
plot(shift_samples,Rdd)
title('Circular Autocorrelation Rdd(tau)');
xlabel('Samples');
ylabel('Rdd');

figure;
fchip=Rc_L1_CA;
samples_per_chip=Fs/fchip;
shift_chips=shift_samples/samples_per_chip;
plot(shift_chips,Rdd)
title('Circular Autocorrelation Rdd(tau)');
xlabel('Chips');
ylabel('Rdd');

% **** QUESTION 4 ****:
% Compute the estimation of the autocorrelation gain defined as the ratio 
% between abs(Rdd(tau=0)) and abs(Rdd(tau>>1chip)). Convert it to dB.

% **** QUESTION 5 ****:
% Compute the estimation of the cross-correlation gain between different satellite's PRNs defined as the ratio 
% between abs(R_(d_sat1,d_sat1)(tau=0)) and abs(R_(d_sat1,d_sat2)(tau=0)). Convert it to dB.

%% Autocorrelation vs. frequency shift
%autocorrelation maximum abs(Rdd(tau=0)) vs. frequency shift

d=s_bb_sat(1,:);
Fd_error_max_hz=2000;
shift_hz=-Fd_error_max_hz:10:Fd_error_max_hz;

% **** QUESTION 6 ****:
% compute the autocorrelation estimation of Rdd(tau=0,f) vs. a frequency
% doppler error simulated as a frequency shift.
% HINT: R_(d,d(f))=(1/K)*(d*d*exp(-jW))'
% where w is expressed in rad/s

for n=xxxxx
%     Rdd_doppler_error(n)=?;
end
pause;

% **** QUESTION 7 ****:
% Measure the 50% autocorrelation gain reduction vs. frequency. Express the
% result in Hz. Define the Doppler shift acquisition grid granularity based on the
% result.
% Increase now the signal duration Tsig from 1ms to several ms. Plot again
% the is sinc-shaped curve. Recompute the 50% gain reduction and the
% Doppler shift grid granularity and discuss the results.

figure;
plot(shift_hz,abs(Rdd_doppler_error))
title('Circular Autocorrelation Rdd(tau=0,f)');
xlabel('Freq error [Hz]');
ylabel('abs(Rdd)');



