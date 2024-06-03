% QPSK modulation
% Symbol rate (Fs) 1 ksps => symbol interval (Ts) 1 ms
% Pulse shaping filter: roll-off factor 0.4, oversampling (RATE) 10
% Carrier frequency : 2 kHz
close all; clear all; clf;
Ts=10^(-3); % symbol interval 1 ms
no_data=200000;
temp1=rand(1,no_data); temp2=rand(1,no_data); %Uniformly distributed random numbers
format long
% Bit data generation

P_dB = 1:17;

for i1 = 1:length(P_dB)

    p = 10^(P_dB(i1)*0.1);

    for i=1:no_data
        if(temp1(i)<0.5) int1(i)=-1;
        else int1(i)=1;
        end
        if(temp2(i)<0.5) int2(i)=-1;
        else int2(i)=1;
        end
    end

    sym1 = int1*sqrt(p/2);
    sym2 = int2*sqrt(p/2);

    
    % Pulse shaping filtering: raised cosine filter
    roll_off=0.4; N_T=5; RATE=10; Fs=1/(Ts/RATE);
    p=rcosdesign(roll_off, N_T, RATE, 'sqrt');
    
    
    
    % Data to the psf
    for i=1:no_data
        int_sym1((i-1)*RATE+1)=sym1(i);
        int_sym2((i-1)*RATE+1)=sym2(i);
        for j=2:RATE
            int_sym1((i-1)*RATE+j)=0;
            int_sym2((i-1)*RATE+j)=0;
        end
    end
    
    % Generation of pulse shaped signal
    sym_p1=conv(p,int_sym1);
    sym_p2=conv(p,int_sym2);
    
    % Demodulated signal
    rx1=sym_p1 + randn(1, length(sym_p1));
    rx2=sym_p2 + randn(1, length(sym_p2));
    
    % Matched filter
    rec_p1=conv(p,rx1);
    rec_p2=conv(p,rx2);
    
    % Switch & Detect
    for i = 1: no_data
        yI(i) = rec_p1((i+9)*10);
        yQ(i) = rec_p2((i+9)*10);

        if yI(i)>0
            yI_detect(i) = 1;
        else
            yI_detect(i) = -1;
        end

        if yQ(i)>0
            yQ_detect(i) = 1;
        else
            yQ_detect(i) = -1;
        end
    end

    SER(i1) = 1- nnz((int1 == yI_detect)&(int2 == yQ_detect))/no_data;
end

semilogy(P_dB, SER)

hold on;

pp = 10.^(P_dB(1:17)*0.1);
er = qfunc(sqrt(pp/2));
semilogy(P_dB(1:17), er);

legend('Simulation', 'Theoretical')


P_dB = [10 15 20 30];
for i1 = 1:length(P_dB)

    p = 10^(P_dB(i1)*0.1);

    for i=1:no_data
        if(temp1(i)<0.5) int1(i)=-1;
        else int1(i)=1;
        end
        if(temp2(i)<0.5) int2(i)=-1;
        else int2(i)=1;
        end
    end

    sym1 = int1*sqrt(p/2);
    sym2 = int2*sqrt(p/2);

    
    % Pulse shaping filtering: raised cosine filter
    roll_off=0.4; N_T=5; RATE=10; Fs=1/(Ts/RATE);
    p=rcosdesign(roll_off, N_T, RATE, 'sqrt');
    
    
    
    % Data to the psf
    for i=1:no_data
        int_sym1((i-1)*RATE+1)=sym1(i);
        int_sym2((i-1)*RATE+1)=sym2(i);
        for j=2:RATE
            int_sym1((i-1)*RATE+j)=0;
            int_sym2((i-1)*RATE+j)=0;
        end
    end
    
    
    % Generation of pulse shaped signal
    sym_p1=conv(p,int_sym1);
    sym_p2=conv(p,int_sym2);
    
    
    
    % Carrier signal generation
    Fc=2000; % Carrier frequency
    N=length(sym_p1);
    t=Ts/RATE*(1:N);
    C1=sqrt(2)*cos(2*pi*Fc*t);
    C2=sqrt(2)*sin(2*pi*Fc*t);
    
    
    % Modulated signal
    tx=sym_p1.*C1 + sym_p2.*C2;
    
    % AWGN
    rx = tx+randn(1, length(tx));

    % Demodulated signal
    rx1=rx.*C1;
    rx2=rx.*C2;
    
    % Matched filter
    rec_p1=conv(p,rx1);
    rec_p2=conv(p,rx2);
    
    % Switch & Detect
    for i = 1: no_data
        yI(i) = rec_p1((i+9)*10);
        yQ(i) = rec_p2((i+9)*10);

    end
    
    figure();
    scatter(yI, yQ, 'Marker', '*');

end
