% BPSK modulation
% Symbol rate (Fs) 1 ksps => symbol interval (Ts) 1 ms
% Pulse shaping filter: roll-off factor 0.4, oversampling (RATE) 10
% Carrier frequency : 2 kHz
close all; clear; clc;

Ts = 10^(-3); % symbol interval 1 ms
no_data = 2000; 
sym = 2*randi([0 1],1, no_data)-1; % Randomly distributed -1, 1 bits

% Data spectrum
F_data = fft(sym(1:2000));
F_data = fftshift(F_data); % Shift zero-frequency component to center of spectrum
a = max(abs(F_data));
F_data = F_data / a; % 최대 크기로 스케일링

figure(1);
N = no_data;
freq = (-N/2:N/2-1) / (Ts * N); % frequency vector
plot(freq(N/2:N), 20*log10(abs(F_data(N/2:N)))), grid % plot magnitude spectrum
xlabel('Frequency [Hz]'); ylabel('Normalized Magnitude [dB]')
text(150, -50, 'Sampling Frequency: 1KHz', 'Color', 'red', 'FontSize', 14)
axis([0 500 -60 0])

% Pulse shaping filtering: raised cosine filter
roll_off = 0.4; 
N_T = 5; 
RATE = 10;

% 기존 코드: p = rcosfir(roll_off, N_T,           RATE,                 Ts,                 'sqrt');
%                                  Filter extent  oversampling rate     input sample rate

p = rcosdesign(roll_off, N_T, RATE, 'sqrt');
t = -5*Ts:Ts/RATE:5*Ts;
figure(2); stem(t*1000, p); grid;
xlabel('Time [ms]'); ylabel('Amplitude'); axis([-5 5 -0.1 0.4])

% Data to the pulse shaping filter
for i = 1:no_data
    int_sym((i-1)*RATE+1) = sym(i);
    for j = 2:RATE
        int_sym((i-1)*RATE+j) = 0;
    end
end

% Generation of pulse shaped signal
sym_p = conv(p, int_sym);
figure(3); plot(sym_p(1100:1300)); grid;
xlabel('Time'); 

% Eye diagram
t_int = 0:Ts/RATE*1000:Ts*(2-1/RATE)*1000; % 20 sample interval
figure(4); hold on;
for i = 100:200
    plot(t_int, sym_p((i-1)*RATE*2+1:i*RATE*2))
end
hold off;
xlabel('Time [ms]'); ylabel('Amplitude'); grid;

% Carrier signal generation
Fc = 2000; % Carrier frequency
N = length(sym_p);
t = Ts/RATE*(1:N);
C = cos(2*pi*Fc*t);

% Baseband spectrum plot
freq = (-N/2:N/2-1) / (Ts/RATE * N);  % frequency vector
tmp_Fre = fft(sym_p(1:N));  % DFT/FFT
Fsym_p = fftshift(tmp_Fre); % shift it for plotting
a = max(abs(Fsym_p));
Fsym_p = Fsym_p / a;
figure(5);
plot(freq(N/2:N), 20*log10(abs(Fsym_p(N/2:N)))), grid 
xlabel('Frequency [Hz]'); ylabel('Normalized Magnitude [dB]')
axis([0 5000 -120 0])

% Modulated signal
tx = sym_p .* C;

% Passband spectrum plot
tmp_Fre = fft(tx(1:N)); % DFT/FFT
Ftx = fftshift(tmp_Fre); % shift it for plotting
a = max(abs(Ftx));
Ftx = Ftx / a;
figure(6);
plot(freq(N/2:N), 20*log10(abs(Ftx(N/2:N)))), grid
xlabel('Frequency [Hz]'); ylabel('Normalized Magnitude [dB]')
axis([0 5000 -120 0])

% Demodulated signal
rx = tx .* C;

% Matched filter 
rec_p = conv(p, rx);
figure(7); plot(rec_p(1100:1300)); grid;
xlabel('Time'); 

% Eye diagram
t_int = 0:Ts/RATE*1000:Ts*2*1000;
figure(8); hold on;
for i = 100:200
    plot(t_int, rec_p((i-1)*RATE*2+1:i*RATE*2+1))
end
hold off;
xlabel('Time [ms]'); ylabel('Amplitude'); grid;