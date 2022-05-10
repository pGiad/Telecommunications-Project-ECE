% PANAGIOTIS GIADIKIAROGLOU
% AM 03119185 -> 1 + 8 + 5 = 14 -> 1 + 4 = 5
% -------------------------------------------------------------------------
% MICHAELANGELO VELALOPOULOS
% AM 03119908 -> 9 + 0 + 8 = 17 -> 1 + 7 = 8

% Group 81

%--------------------------------------------------------------------------

% Code with Panagiotis Giadikiaroglou parameters


% 1o Erwthma

A = 4;
fm = 5000;
T = 4*(1/fm);
fs = 1000000;
t = 0:1/fs:T-1/fs;
x = A * (sawtooth(2*pi*fm*t, 1/2));
figure (1)
plot(t,x)
grid on
title('Triangular Pulse')
xlabel('Time [s]')
ylabel('Amplitude [V]')

% 1.a.i
fs1 = 30 * fm;
dt1 = 1 / fs1;
t1 = 0 : dt1 : T-1/fs;
x1 = A * (sawtooth(2*pi*fm*t1, 1/2));
figure (2)
stem(t1, x1, 'rx-')
grid on 
title('Sampled Pulse (fs1 = 30 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

% 1.a.ii
fs2 = 50 * fm;
dt2 = 1 / fs2;
t2 = 0 : dt2 : T-1/fs;
x2 = A * (sawtooth(2*pi*fm*t2, 1/2));
figure (3)
stem(t2, x2, 'bo-')
grid on 
title('Sampled Pulse (fs2 = 50 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

% 1.a.iii Ypoerwthma
figure (4)
hold on
stem(t1, x1, 'rx-')
stem(t2, x2, 'bo-')
title('Both Sampled Pulses')
xlabel('Time [s]')
ylabel('Amplitude [V]')
legend('fs1 = 30*fm','fs2 = 50*fm');
grid on
hold off


% 1.b
% DFT of the triangular pulse 
Y = fft(x);
N = length(Y);
f_y = 0:fs/N:fs/2-fs/N;
figure (5)
plot(f_y, abs(Y(1:N/2)));
grid on
title('Phase spectrum of the triangular pulse')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 200000])


fs_b = 4 * fm;
dt = 1 / fs_b;
t_b = 0 : dt : T-1/fs;
x_b = A * (sawtooth(2*pi*fm*t_b, 1/2));
figure (6)
stem(t_b, x_b, 'bo-')
grid on 
title('Sampled Pulse (fs = 4 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')


% 1.c
z = sin(2*pi*fm*t);
figure (7)
plot(t, z);
grid on
title('z(t) = sin(2πfmt)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

% 1.c.i
% 1.c.i.a
z1 = sin(2*pi*fm*t1);
figure (8)
stem(t1, z1, 'rx-')
grid on 
title('Sampled z Signal (fs1 = 30 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

z2 = sin(2*pi*fm*t2);
figure (9)
stem(t2, z2, 'bo-')
grid on 
title('Sampled z Signal (fs2 = 50 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

figure (10)
hold on
stem(t1, z1, 'rx-')
stem(t2, z2, 'bo-')
title('Both Sampled z Signals')
xlabel('Time [s]')
ylabel('Amplitude [V]')
legend('fs1 = 30*fm','fs2 = 50*fm');
grid on
hold off

% 1.c.i.b
z_b = sin(2*pi*fm*t_b);
figure (11)
stem(t_b, z_b, 'bo-')
grid on 
title('Sampled z Signal (fs = 4 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

% 1.c.ii
fm_total = gcd(fm, fm+1000);
T_1 = (1/fm_total);
t_1 = 0:1/fs:T_1-1/fs;
q = sin(2*pi*fm*t_1) + sin(2*pi*(fm + 1000)*t_1);
figure (12)
plot(t_1, q)
grid on
title('q(t) = z(t) + sin(2π(fm+1000)t)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

% 1.c.ii.a
t1_c = 0 : dt1 : T_1-1/fs;
q1 = sin(2*pi*fm*t1_c) + sin(2*pi*(fm + 1000)*t1_c);
figure (13)
stem(t1_c, q1, 'rx-')
grid on 
title('Sampled q Signal (fs1 = 30 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

t2_c = 0 : dt2 : T_1-1/fs;
q2 = sin(2*pi*fm*t2_c) + sin(2*pi*(fm + 1000)*t2_c);
figure (14)
stem(t2_c, q2, 'bo-')
grid on 
title('Sampled q Signal (fs2 = 50 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

figure (15)
hold on
stem(t1_c, q1, 'rx-')
stem(t2_c, q2, 'bo-')
title('Both Sampled q Signals')
xlabel('Time [s]')
ylabel('Amplitude [V]')
legend('fs1 = 30*fm','fs2 = 50*fm');
grid on
hold off

% 1.c.ii.b
t_c_ii_b = 0 : dt : T_1-1/fs;
q_b = sin(2*pi*fm*t_c_ii_b) + sin(2*pi*(fm + 1000)*t_c_ii_b);
figure (16)
stem(t_c_ii_b, q_b, 'bo-')
grid on 
title('Sampled q Signal (fs = 4 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')


%--------------------------------------------------------------------------

% 2o Erwthma

Max_Value_Dec = 2^5; %levels
dec_sequence = (0:Max_Value_Dec-1);
gray_sequence = qammod(dec_sequence, Max_Value_Dec, 'gray');

[gray_sec,map_gray] = bin2gray(dec_sequence,'qam',Max_Value_Dec);
gray_code = dec2bin(gray_sec);


% 2.a
step_size = (max(x1) - min(x1))/Max_Value_Dec; %step size
index_2a = round((x1 - min(x1))/step_size) ; %index
quants = min(x1) + index_2a*step_size;
figure (17)
stem(t1,quants)
grid on
title('Quantized Signal with 5 bits')
xlabel('Time [s]')
ylabel('Gray Code')
yticks([-4 -3.75 -3.5 -3.25 -3 -2.75 -2.5 -2.25 -2 -1.75 -1.5 -1.25 -1 -0.75 -0.5 -0.25 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4])
yticklabels(gray_code);
% ta 4 bits einai arketa gia na anaparastathei olo to sima upo
% deigmatolipsia kai etsi ta 5 bits den exoun kapoia diafora sta epipeda


% 2.b.i
error = x1 - quants;
error_first10 = error(1,1:10);
std_first10 = std(error_first10);
% Result = 0.068

% 2.b.ii
error_first20 = error(1,1:20);
std_first20 = std(error_first20);
% Result = 0.0749

% 2.b.iii
numerator = 3*2^(2*5);
f1 = max(x1(1:10))/rms(x1(1:10));
snr10 = 10*log10(numerator./(f1^2));
% snr10 = 43.7203
f2 = max(x1(1:20))/rms(x1(1:20));
snr20 = 10*log10(numerator./(f2^2));
% snr20 = 30.8463
f_x1 = max(x1)/rms(x1);
snr_x1 = 10*log10(numerator./(f_x1^2));
% snr_x1 = 30.1414


% 2.c
bit_stream = char(zeros(30,5));
for i = 1:30
    if quants(i) < 0
        index_gc = (quants(i) + 4) / 0.25 + 1;
        bit_stream(i,:) = gray_code(index_gc,:);
    else 
        index_gc = (quants(i) + 4) / 0.25;
        bit_stream(i,:) = gray_code(index_gc,:);
    end
end

x_axis_bit_stream = [];
y_axis_bit_stream = [];
index_bs = 1;
index_xy = 1;
for i = 1:30
    for j = 1:5
        x_axis_bit_stream(index_xy:index_xy+3) =  [2*(index_bs-1) 2*(index_bs-1+0.5) 2*(index_bs-1+0.5) 2*(index_bs)];
        if(bit_stream(i,j) == '0')
            y_axis_bit_stream(index_xy:index_xy+3) = [-5 -5 0 0];
            index_xy = index_xy + 4;
        else 
            y_axis_bit_stream(index_xy:index_xy+3) = [5 5 0 0];
            index_xy = index_xy + 4;
        end
        index_bs = index_bs + 1;
    end
end

figure (18)
plot(x_axis_bit_stream, y_axis_bit_stream), axis([0, 300, -6, 6]);
grid on
title('Polar RZ');
xlabel('Time [ms]')
ylabel('Amplitude [V]')


%--------------------------------------------------------------------------

% 3o Erwthma

random_bits = randi([0,1],[36,1]);

bits_length = length(random_bits);
Tb = 0.25;
fc = 2;
nb = 100;

random_bits_string = [];
for i = 1:bits_length
    if random_bits(i) == 1
        bit_signal = ones(1, nb);
    else 
        bit_signal = zeros(1, nb);
    end
    random_bits_string = [random_bits_string bit_signal];
end

t_random_bits = Tb/nb:Tb/nb:bits_length*Tb;
figure (19)
plot(t_random_bits, random_bits_string);
grid on
xlabel('Time [s]')
ylim([-0.2 1.2])
ylabel('Amplitude [V]')
title('Digital Signal of Random Bit String')  

% 3.a

% BPSK Modulation
t_bpsk = Tb/nb:Tb/nb:Tb;
mod_bpsk = [];
for i = 1:bits_length
    if random_bits(i) == 1
        bpsk_bits = cos(2*pi*fc*t_bpsk);
    else
        bpsk_bits = cos(2*pi*fc*t_bpsk + pi);
    end
    mod_bpsk = [mod_bpsk bpsk_bits];
end

figure (20)
plot(t_random_bits, mod_bpsk);
grid on
xlabel('Time [s]')
ylabel('Amplitude [V]')
title('BPSK Modulated Signal')

% QPSK Modulation
t_qpsk = Tb/nb:Tb/nb:2*Tb;
mod_qpsk = [];
for i = 1:2:bits_length-1
    if random_bits(i) == 0
        if random_bits(i+1) == 0
            qpsk_bits = cos(2*pi*fc*t_qpsk);
        else
            qpsk_bits = sin(2*pi*fc*t_qpsk);
        end
    else 
        if random_bits(i+1) == 0
            qpsk_bits = -cos(2*pi*fc*t_qpsk);            
        else
            qpsk_bits = -sin(2*pi*fc*t_qpsk);
        end
    end
    mod_qpsk = [mod_qpsk qpsk_bits];
end

figure (21)
plot(t_random_bits, mod_qpsk);
grid on
xlabel('Time [s]')
ylabel('Amplitude [V]')
title('QPSK Modulated Signal')

% 8PSK Modulation
t_8psk = Tb/nb:Tb/nb:3*Tb;
mod_8psk = [];
for i = 1:3:bits_length-2
    if random_bits(i) == 0
        if random_bits(i+1) == 0
            if random_bits(i+2) == 0
                bits_8psk = cos(2*pi*fc*t_8psk);
            else
                bits_8psk = cos(2*pi*fc*t_8psk + pi/4);
            end
        else
            if random_bits(i+2) == 0
                bits_8psk = cos(2*pi*fc*t_8psk + 3*pi/4);
            else   
                bits_8psk = cos(2*pi*fc*t_8psk + pi/2);
            end
        end
    else 
        if random_bits(i+1) == 0
            if random_bits(i+2) == 0
                bits_8psk = cos(2*pi*fc*t_8psk + 7*pi/4);
            else
                bits_8psk = cos(2*pi*fc*t_8psk + 3*pi/2);
            end
        else
            if random_bits(i+2) == 0
                bits_8psk = cos(2*pi*fc*t_8psk + pi);
            else
                bits_8psk = cos(2*pi*fc*t_8psk + 5*pi/4);
            end
        end
    end
    mod_8psk = [mod_8psk bits_8psk];
end

figure (22)
plot(t_random_bits, mod_8psk);
grid on
xlabel('Time [s]')
ylabel('Amplitude [V]')
title('8-PSK Modulated Signal')

% 3.b
A_3b = 5;
bpam_random_bits = 2*A_3b*(random_bits_string - 0.5);
figure (23)
plot(t_random_bits, bpam_random_bits);
grid on
xlabel('Time [s]')
ylabel('Amplitude [V]')
title('B-PAM Signal')
ylim([-6 6])

% 3.c
scatterplot(0.5*bpam_random_bits, 1, 0, 'y*');
grid on
title('B-Pam Constellation Diagram')

% 3.d
awgn_3d_5db = awgn(bpam_random_bits, 5);
figure (25)
plot(t_random_bits, awgn_3d_5db)
grid on
xlabel('Time [s]')
ylabel('Amplitude [V]')
title('AWGN for E_b/N_0 = 5dB')

awgn_3d_15db = awgn(bpam_random_bits, 15);
figure (26)
plot(t_random_bits, awgn_3d_15db)
grid on
xlabel('Time [s]')
ylabel('Amplitude [V]')
title('AWGN for E_b/N_0 = 15dB')

% 3.e
bpam_mod_3e = pammod(random_bits, 2);

awgn_bpam_mod_5dB = awgn(bpam_mod_3e, 5+3);
scatterplot((A_3b/2)*awgn_bpam_mod_5dB, 1, 0, 'y*');
grid on
title('B-PAM with E_b/N_0 = 5dB')

awgn_bpam_mod_15dB = awgn(bpam_mod_3e, 15+3);
scatterplot((A_3b/2)*awgn_bpam_mod_15dB, 1, 0, 'y*');
grid on
title('B-PAM with E_b/N_0 = 15dB')


% 3.st
bit_size = 2;
EbN0 = 0:1:15;
[BER] = berawgn(EbN0, 'pam', bit_size);

n_bits_3st = 1000000;
k_3st = log2(bit_size);
snr_dB = EbN0+3+log10(k_3st);
y_noisy = zeros(n_bits_3st,length(snr_dB));
z_snr = zeros(n_bits_3st,length(snr_dB));
errVec = zeros(3,length(EbN0));
errcalc = comm.ErrorRate;
random_bits_3st = randi([0,1],[1000000,1]);
y_random_bits_pam = pammod(random_bits_3st, bit_size);

for jj = 1:length(snr_dB)
    reset(errcalc)
    y_noisy(:,jj) = awgn(real(y_random_bits_pam), snr_dB(jj),'measured');
    z_snr(:,jj) = pamdemod(complex(y_noisy(:,jj)),bit_size);
    errVec(:,jj) = errcalc(random_bits_3st, z_snr(:,jj));
end

figure (29)
semilogy(EbN0, BER, 'r');
hold on
semilogy(EbN0,errVec(1,:),'b.');
grid on
legend('Theoretical BER', 'Experimental BER')
title('Theoretical-Experimental BER for B-PAM')
xlabel('E_b/N_0')
ylabel('Bit Error Rate')
hold off


%--------------------------------------------------------------------------

% 4o Erwthma

% 4.a
qpsk_4 = comm.QPSKModulator(BitInput=true);
qpsk_bits_4 = A_3b*qpsk_4(random_bits);
scatterplot(0.5*qpsk_bits_4, 1, 0, 'y*');
grid on
title('QPSK Modulated Signal')
text(0.5*A_3b/sqrt(2)-0.12, 0.5*A_3b/sqrt(2)-0.22, '00', 'Color', [1 1 1]);
text(-0.5*A_3b/sqrt(2)-0.12, 0.5*A_3b/sqrt(2)-0.22, '01', 'Color', [1 1 1]);
text(-0.5*A_3b/sqrt(2)-0.12, -0.5*A_3b/sqrt(2)+0.22, '11', 'Color', [1 1 1]);
text(0.5*A_3b/sqrt(2)-0.12, -0.5*A_3b/sqrt(2)+0.22, '10', 'Color', [1 1 1]);


% 4.b
awgn_4b_5db = awgn(qpsk_bits_4, 5);
scatterplot(0.5*awgn_4b_5db, 1, 0, 'y*');
grid on
title('QPSK with E_b/N_0 = 5dB')

awgn_4b_15db = awgn(qpsk_bits_4, 15);
scatterplot(0.5*awgn_4b_15db, 1, 0, 'y*');
grid on
title('QPSK with E_b/N_0 = 15dB')


% 4.c
bit_size_qpsk = 4;
EbN0 = 0:1:15;
[BER_4c] = berawgn(EbN0,'psk',bit_size_qpsk ,'diff');
figure(33)
semilogy(EbN0,BER_4c,'b')
grid on
title('Theoretical-Experimental BER for QPSK and BPSK')
xlabel('E_b/N_0')
ylabel('Bit Error Rate')

bit_size_bpsk = 2;
n_bits_bpsk = 1000000;
k_bits_bpsk = log2(bit_size_bpsk);
snr_random_bits_bpsk = EbN0+10*log10(k_bits_bpsk);
y_noisy_bpsk = zeros(n_bits_bpsk,length(snr_random_bits_bpsk));
z_snr_bpsk = zeros(n_bits_bpsk,length(snr_random_bits_bpsk));
errVec_bpsk = zeros(3,length(EbN0));

errcalc = comm.ErrorRate;

x_random_bits_bpsk = randi([0 1],1000000,1);
y_random_bits_bpsk = pskmod(x_random_bits_bpsk, bit_size_bpsk);
for jj = 1:length(snr_random_bits_bpsk)
    reset(errcalc)
    y_noisy_bpsk(:,jj) = awgn(y_random_bits_bpsk,snr_random_bits_bpsk(jj),'measured');
    z_snr_bpsk(:,jj) = pskdemod(complex(y_noisy_bpsk(:,jj)), bit_size_bpsk);
    errVec_bpsk(:,jj) = errcalc(x_random_bits_bpsk,z_snr_bpsk(:,jj));
end

hold on
semilogy(EbN0,errVec_bpsk(1,:),'g.');

n_bits_qpsk = 1000000;
k_bits_qpsk = log2(bit_size_qpsk);
snr_random_bits_qpsk = EbN0+10*log10(k_bits_qpsk);
y_noisy_qpsk = zeros(n_bits_qpsk,length(snr_random_bits_qpsk));
z_snr_qpsk = zeros(n_bits_qpsk,length(snr_random_bits_qpsk));
errVec_qpsk = zeros(3,length(EbN0));

errcalc = comm.ErrorRate;

x_random_bits_qpsk = randi([0 1],1000000,1);
y_random_bits_qpsk = pskmod(x_random_bits_qpsk, bit_size_qpsk);

for jj = 1:length(snr_random_bits_qpsk)
    reset(errcalc)
    y_noisy_qpsk(:,jj) = awgn(y_random_bits_qpsk,snr_random_bits_qpsk(jj),'measured');
    z_snr_qpsk(:,jj) = pskdemod(complex(y_noisy_qpsk(:,jj)), bit_size_qpsk);
    errVec_qpsk(:,jj) = errcalc(x_random_bits_qpsk,z_snr_qpsk(:,jj));
end

semilogy(EbN0,errVec_qpsk(1,:),'r.');
legend('Theoretical BER','BPSK Experimental BER','QPSK Experimental BER');
hold off


% 4.d.i
text_file = fopen('rice_odd.txt', 'r');
text_content = fscanf(text_file, '%c');
text_bits = dec2bin(text_content);

% 4.d.ii
text_bits_dec = bin2dec(text_bits);
text_bits_quant = 8;
text_levels = 2^text_bits_quant;
text_step_size = (max(text_bits_dec) - min(text_bits_dec))/text_levels;
text_index = round((text_bits_dec - min(text_bits_dec))/text_step_size);
quantized_text = min(text_bits_dec) + text_index*text_step_size;
figure (35);
n_text = 1:1:485;
stem(n_text, quantized_text, 'rx-');
grid on
xlabel('Characters')
ylabel('Decimal Representation of Ascii Code')
title('Quantized Text Signal')

% 4.d.iii
quantized_text_bits = dec2bin(round(quantized_text)) - '0';
quantized_bit_stream = reshape(quantized_text_bits.', [], 1);
quantized_bit_stream(end+1) = 0;

text_qpsk_mod = qpsk_4(quantized_bit_stream);
scatterplot(text_qpsk_mod, 1, 0, 'y*');
grid on
title('Text QPSK Constellation Diagram')

% 4.d.iv
awgn_text_qpsk_mod_5dB = awgn(text_qpsk_mod, 5);
awgn_text_qpsk_mod_15dB = awgn(text_qpsk_mod, 15);

% 4.d.v
text_qpsk_demod = comm.QPSKDemodulator(BitOutput=true); 
text_qpsk_demod_5dB = text_qpsk_demod(awgn_text_qpsk_mod_5dB);
text_qpsk_demod_15dB = text_qpsk_demod(awgn_text_qpsk_mod_15dB);
scatterplot(awgn_text_qpsk_mod_5dB);
grid on
title('Text QPSK with E_s/N_0 = 5dB')
scatterplot(awgn_text_qpsk_mod_15dB);
grid on
title('Text QPSK with E_s/N_0 = 15dB')

% 4.d.vi
[num_qpsk_demod_5dB, error_qpsk_demod_5dB] = symerr(quantized_bit_stream, text_qpsk_demod_5dB);
[num_qpsk_demod_15dB, error_qpsk_demod_15dB] = symerr(quantized_bit_stream, text_qpsk_demod_15dB);
ber_theoretical_5dB = qfunc(sqrt(10^(5/10)));
% ber_theoretical_5dB = 0.0377
ber_theoretical_15dB = qfunc(sqrt(10^(15/10)));
% ber_theoretical_15dB = 0

% 4.d.vii
text_demod_5dB  = char(bin2dec(reshape(char('0' + text_qpsk_demod_5dB(1:3395)),7,[]).'))';
text_demod_15dB  = char(bin2dec(reshape(char('0' + text_qpsk_demod_15dB(1:3395)),7,[]).'))';
text_5dB = fopen('text_5db_03119185.txt', 'w');
fprintf(text_5dB, '%c', text_demod_5dB);
text_15dB = fopen('text_15db_03119185.txt', 'w');
fprintf(text_15dB, '%c', text_demod_15dB);


%--------------------------------------------------------------------------

% 5o Erwthma

% 5.a
[sound_wave, fs_sound] = audioread('soundfile1_lab2.wav'); 
dt_sound = 1/fs_sound;
t_sound = 0:dt_sound:(length(sound_wave)*dt_sound)-dt_sound;
figure (39);
plot(t_sound, sound_wave);
grid on
xlabel('Time [s]')
ylabel('Amplitude')
title('Sound Signal')


% 5.b
sound_bits = 8;
sound_levels = 2^sound_bits;
sound_step_size = (max(sound_wave) - min(sound_wave))/sound_levels;
sound_index = round((sound_wave - min(sound_wave))/sound_step_size);
quantized_sound = min(sound_wave) + sound_index * sound_step_size;
figure (40);
stem(t_sound,quantized_sound, 'gx-');
grid on
xlabel('Time [s]')
ylabel('Amplitude')
title('Quantized Sound Signal at 8 bits')


% 5.c
quantized_sound_mult = quantized_sound * 1000;
quantized_sound_bits = dec2bin(quantized_sound_mult,16) - '0';
quantized_sound_bit_stream = reshape(quantized_sound_bits.',[],1);
qpsk_mod_5 = comm.QPSKModulator('BitInput', true);
sound_qpsk_mod = qpsk_mod_5(quantized_sound_bit_stream);
scatterplot(sound_qpsk_mod, 1, 0, 'y*');
grid on
title('Sound QPSK Constellation Diagram')


% 5.d
awgn_sound_qpsk_mod_4dB = awgn(sound_qpsk_mod, 4);
awgn_sound_qpsk_mod_14dB = awgn(sound_qpsk_mod, 14);


% 5.e
qpsk_demod_5 = comm.QPSKDemodulator('BitOutput',true); 
sound_qpsk_demod_4dB = qpsk_demod_5(awgn_sound_qpsk_mod_4dB);
sound_qpsk_demod_14dB = qpsk_demod_5(awgn_sound_qpsk_mod_14dB);
scatterplot(awgn_sound_qpsk_mod_4dB);
grid on
title('Sound QPSK with E_s/N_0 = 4dB')
scatterplot(awgn_sound_qpsk_mod_14dB);
grid on
title('Sound QPSK with E_s/N_0 = 14dB')


% 5.st
[num_qpsk_demod_4dB, error_qpsk_demod_4db] = symerr(quantized_sound_bit_stream, sound_qpsk_demod_4dB);
[num_qpsk_demod_14dB, error_qpsk_demod_14db] = symerr(quantized_sound_bit_stream, sound_qpsk_demod_14dB);
ber_theoretical_4dB = qfunc(sqrt(10^(4/10)));
% ber_theoretical_4dB = 0.0565
ber_theoretical_14dB = qfunc(sqrt(10^(14/10)));
% ber_theoretical_14dB = 0


% 5.z
sound_reshaped_14dB = char('0'+ (reshape(sound_qpsk_demod_14dB, 16, []).'));
sound_demod_dec_14dB = typecast(uint16(bin2dec(sound_reshaped_14dB)),'int16');
sound_demod_14dB = double(abs(sound_demod_dec_14dB))./1000 ;

sound_reshaped_4dB = char('0'+ (reshape(sound_qpsk_demod_4dB, 16, []).'));
sound_demod_dec_4dB = typecast(uint16(bin2dec(sound_reshaped_4dB)),'int16');
sound_demod_4dB = double(abs(sound_demod_dec_4dB))./1000 ;

sound_demod_14dB_norm = sound_demod_14dB / max(sound_demod_14dB);
audiowrite("sound_14dB_03119185.wav", sound_demod_14dB_norm , 44100)
sound_demod_4dB_norm = sound_demod_4dB / max(sound_demod_4dB);
audiowrite("sound_4dB_03119185.wav", sound_demod_4dB_norm, 44100)


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


% Code with Michaelangelo Velalopoulos parameters

%{
% 1o Erwthma

A = 4;
fm = 8000;
T = 4*(1/fm);
fs = 1000000;
t = 0:1/fs:T-1/fs;
x = A * (sawtooth(2*pi*fm*t, 1/2));
figure (1)
plot(t,x)
grid on
title('Triangular Pulse')
xlabel('Time [s]')
ylabel('Amplitude [V]')

% 1.a.i
fs1 = 30 * fm;
dt1 = 1 / fs1;
t1 = 0 : dt1 : T-1/fs;
x1 = A * (sawtooth(2*pi*fm*t1, 1/2));
figure (2)
stem(t1, x1, 'rx-')
grid on 
title('Sampled Pulse (fs1 = 30 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

% 1.a.ii
fs2 = 50 * fm;
dt2 = 1 / fs2;
t2 = 0 : dt2 : T-1/fs;
x2 = A * (sawtooth(2*pi*fm*t2, 1/2));
figure (3)
stem(t2, x2, 'bo-')
grid on 
title('Sampled Pulse (fs2 = 50 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

% 1.a.iii Ypoerwthma
figure (4)
hold on
stem(t1, x1, 'rx-')
stem(t2, x2, 'bo-')
title('Both Sampled Pulses')
xlabel('Time [s]')
ylabel('Amplitude [V]')
legend('fs1 = 30*fm','fs2 = 50*fm');
grid on
hold off


% 1.b
% DFT of the triangular pulse 
Y = fft(x);
N=length(Y);
f_y = 0:fs/N:fs/2-fs/N;
figure (5)
plot(f_y, abs(Y(1:N/2)));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Phase spectrum of the triangular pulse');

fs_b = 4 * fm;
dt = 1 / fs_b;
t_b = 0 : dt : T-1/fs;
x_b = A * (sawtooth(2*pi*fm*t_b, 1/2));
figure (6)
stem(t_b, x_b, 'bo-')
grid on 
title('Sampled Pulse (fs = 4 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')


% 1.c
z = sin(2*pi*fm*t);
figure (7)
plot(t, z);
grid on
title('z(t) = sin(2πfmt)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

% 1.c.i.a
z1 = sin(2*pi*fm*t1);
figure (8)
stem(t1, z1, 'rx-')
grid on 
title('Sampled z Signal (fs1 = 30 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

z2 = sin(2*pi*fm*t2);
figure (9)
stem(t2, z2, 'bo-')
grid on 
title('Sampled z Signal (fs2 = 50 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

figure (10)
hold on
stem(t1, z1, 'rx-')
stem(t2, z2, 'bo-')
title('Both Sampled z Signals')
xlabel('Time [s]')
ylabel('Amplitude [V]')
legend('fs1 = 30*fm','fs2 = 50*fm');
grid on
hold off

% 1.c.i.b
z_b = sin(2*pi*fm*t_b);
figure (11)
stem(t_b, z_b, 'bo-')
grid on 
title('Sampled z Signal (fs = 4 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

% 1.c.ii
fm_total = gcd(fm, fm+1000);
T_1 = (1/fm_total);
t_1 = 0:1/fs:T_1-1/fs;
q = sin(2*pi*fm*t_1) + sin(2*pi*(fm + 1000)*t_1);
figure (12)
plot(t_1, q)
grid on
title('q(t) = z(t) + sin(2π(fm+1000)t)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

% 1.c.ii.a
t1_c = 0 : dt1 : T_1-1/fs;
q1 = sin(2*pi*fm*t1_c) + sin(2*pi*(fm + 1000)*t1_c);
figure (13)
stem(t1_c, q1, 'rx-')
grid on 
title('Sampled q Signal (fs1 = 30 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

t2_c = 0 : dt2 : T_1-1/fs;
q2 = sin(2*pi*fm*t2_c) + sin(2*pi*(fm + 1000)*t2_c);
figure (14)
stem(t2_c, q2, 'bo-')
grid on 
title('Sampled q Signal (fs2 = 50 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

figure (15)
hold on
stem(t1_c, q1, 'rx-')
stem(t2_c, q2, 'bo-')
title('Both Sampled q Signals')
xlabel('Time [s]')
ylabel('Amplitude [V]')
legend('fs1 = 30*fm','fs2 = 50*fm');
grid on
hold off

% 1.c.ii.b
t_c_ii_b = 0 : dt : T_1-1/fs;
q_b = sin(2*pi*fm*t_c_ii_b) + sin(2*pi*(fm + 1000)*t_c_ii_b);
figure (16)
stem(t_c_ii_b, q_b, 'bo-')
grid on 
title('Sampled q Signal (fs = 4 * fm)')
xlabel('Time [s]')
ylabel('Amplitude [V]')

%--------------------------------------------------------------------------

% 2o Erwthma

Max_Value_Dec = 32;
dec_sequence = (0:Max_Value_Dec-1);
gray_sequence = qammod(dec_sequence, Max_Value_Dec, 'gray');
[gray_sec,map_gray] = bin2gray(dec_sequence,'qam',Max_Value_Dec);
gray_code = dec2bin(gray_sec);


% 2.a
step_size = (max(x1) - min(x1))/Max_Value_Dec; %step size
index_2a = round((x1 - min(x1))/step_size) ; %index
quants = min(x1) + index_2a*step_size;
figure (17)
stem(t1,quants)
grid on
title('Quantized Signal with 5 bits')
xlabel('Time [s]')
ylabel('Gray Code')
yticks([-4 -3.75 -3.5 -3.25 -3 -2.75 -2.5 -2.25 -2 -1.75 -1.5 -1.25 -1 -0.75 -0.5 -0.25 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4])
yticklabels(gray_code);
% ta 4 bits einai arketa gia na anaparastathei olo to sima upo
% deigmatolipsia kai etsi ta 5 bits den exoun kapoia diafora sta epipeda


% 2.b.i
error = x1 - quants;
error_first10 = error(1,1:10);
std_first10 = std(error_first10);
% Result = 0.0046

% 2.b.ii
error_first20 = error(1,1:20);
var_first20 = std(error_first20);
% Result = 0.0056

% 2.b.iii
numerator = 3*2^(2*5);
f1 = max(x1(1:10))/rms(x1(1:10));
snr10 = 10*log10(numerator./(f1^2));
% snr10 = 43.7203
f2 = max(x1(1:20))/rms(x1(1:20));
snr20 = 10*log10(numerator./(f2^2));
% snr20 = 30.8463
f_x1 = max(x1)/rms(x1);
snr_x1 = 10*log10(numerator./(f_x1^2));
%snr_x1 = 30.1414


% 2.c
bit_stream = char(zeros(30,5));
for i = 1:30
    if quants(i) < 0
        index_gc = ( quants(i) + 4 ) / 0.25 + 1;
        bit_stream(i,:) = gray_code(index_gc,:);
    else 
        index_gc = ( quants(i) + 4 ) / 0.25;
        bit_stream(i,:) = gray_code(index_gc,:);
    end
end

x_axis_bit_stream = [];
y_axis_bit_stream = [];
index_bs = 1;
for i = 1:30
    for j = 1:5
        x_axis_bit_stream = [x_axis_bit_stream 2*(index_bs-1) 2*(index_bs-1+0.5) 2*(index_bs-1+0.5) 2*(index_bs)];
        if(bit_stream(i,j) == '0')
            y_axis_bit_stream = [y_axis_bit_stream -8 -8 0 0];
        else 
            y_axis_bit_stream = [y_axis_bit_stream 8 8 0 0];
        end
        index_bs = index_bs + 1;
    end
end

figure (18)
plot(x_axis_bit_stream, y_axis_bit_stream), axis([0, 300, -9, 9]);
title('Polar RZ');
xlabel('Time [ms]')
ylabel('Amplitude [V]')

%--------------------------------------------------------------------------


%3o Erotima

random_bits = randi([0,1],[36,1]);
bs_len = length(random_bits);
Tb = 0.25;
fc = 1;
nb = 100;

rand_bs_digital_string2 = [];
for i = 1:bs_len
    if random_bits(i) == 1
        bit_signal = ones(1, nb);
    else 
        bit_signal = zeros(1, nb);
    end
    rand_bs_digital_string2 = [rand_bs_digital_string2 bit_signal];
end

t_bitstream = Tb/nb:Tb/nb:bs_len*Tb;
figure (19)
grid on
plot(t_bitstream, rand_bs_digital_string2, 'LineWidth', 2.5);
xlabel('Time[s]')
ylim([-0.5 1.5])
ylabel('Amplitude[V]')
title('Digital Signal of Random Bit String')  


% 3.a

% BPSK Modulation
t_1bitperiod = Tb/nb:Tb/nb:Tb;
BPSK_mod = [];
for i = 1:bs_len
    if random_bits(i) == 1
        bit_modulation_BPSK = cos(2*pi*fc*t_1bitperiod);
    else
        bit_modulation_BPSK = cos(2*pi*fc*t_1bitperiod + pi);
    end
    BPSK_mod = [BPSK_mod bit_modulation_BPSK];
end

figure (19) 
grid on
plot(t_bitstream, BPSK_mod, 'LineWidth', 1.5);
xlabel('Time[s]')
ylabel('Amplitude[V]')
title('BPSK Modulated Signal')

% QPSK Modulation
t_2bitperiod = Tb/nb:Tb/nb:2*Tb;
QPSK_mod = [];
for i = 1:2:bs_len-1
    if random_bits(i) == 0
        if random_bits(i+1) == 0
            bit_modulation_QPSK = cos(2*pi*fc*t_2bitperiod);
        else
            bit_modulation_QPSK = sin(2*pi*fc*t_2bitperiod);
        end
    else 
        if random_bits(i+1) == 0
            bit_modulation_QPSK = -cos(2*pi*fc*t_2bitperiod);            
        else
            bit_modulation_QPSK = -sin(2*pi*fc*t_2bitperiod);
        end
    end
    QPSK_mod = [QPSK_mod bit_modulation_QPSK];
end

figure (20)
grid on
plot(t_bitstream, QPSK_mod, 'LineWidth', 1.5);
xlabel('Time[s]')
ylabel('Amplitude[V]')
title('QPSK Modulated Signal')

% 8-PSK Modulation
t_3bitperiod = Tb/nb:Tb/nb:3*Tb;
PSK8_mod = [];
for i = 1:3:bs_len-2
    if random_bits(i) == 0
        if random_bits(i+1) == 0
            if random_bits(i+2) == 0
                bit_modulation_8PSK = cos(2*pi*fc*t_3bitperiod);
            else
                bit_modulation_8PSK = cos(2*pi*fc*t_3bitperiod + pi/4);
            end
        else
            if random_bits(i+2) == 0
                bit_modulation_8PSK = cos(2*pi*fc*t_3bitperiod + 3*pi/4);
            else   
                bit_modulation_8PSK = cos(2*pi*fc*t_3bitperiod + pi/2);
            end
        end
    else 
        if random_bits(i+1) == 0
            if random_bits(i+2) == 0
                bit_modulation_8PSK = cos(2*pi*fc*t_3bitperiod + 7*pi/4);
            else
                bit_modulation_8PSK = cos(2*pi*fc*t_3bitperiod + 3*pi/2);
            end
        else
            if random_bits(i+2) == 0
                bit_modulation_8PSK = cos(2*pi*fc*t_3bitperiod + pi);
            else
                bit_modulation_8PSK = cos(2*pi*fc*t_3bitperiod + 5*pi/4);
            end
        end
    end
    PSK8_mod = [PSK8_mod bit_modulation_8PSK];
end

figure(22)
grid on
plot(t_bitstream, PSK8_mod, 'LineWidth', 1.5);
xlabel('Time[s]')
ylabel('Amplitude[V]')
title('8-PSK Modulated Signal')


% 3.b
A = 8;
B_PAM_rand_bs = 2*A*(rand_bs_digital_string2 - 0.5);
figure(23)
plot(t_bitstream, B_PAM_rand_bs, 'LineWidth', 2.5);
grid on
xlabel('Time[s]')
ylim([-9 9])
ylabel('Amplitude[V]')
title('B-PAM of Random Bits')  


% 3.c
scatterplot(0.5*B_PAM_rand_bs, 1, 0, 'y*');
grid on
title('B-Pam Constellation Diagram')


% 3.d
awgn_3d_5db = awgn(B_PAM_rand_bs, 5);
figure (25)
plot(t_bitstream, awgn_3d_5db)
grid on
xlabel('Time [s]')
ylabel('Amplitude [V]')
title('AWGN for E_b/N_0 = 5db')

awgn_3d_15db = awgn(B_PAM_rand_bs, 15);
figure (26)
plot(t_bitstream, awgn_3d_15db)
grid on
xlabel('Time [s]')
ylabel('Amplitude [V]')
title('AWGN for E_b/N_0 = 15db')


% 3.e
bpam_mod_3e = pammod(random_bits, 2);
awgn_bpam_mod_5dB = awgn(bpam_mod_3e, 5+3);
scatterplot((A/2)*awgn_bpam_mod_5dB, 1, 0, 'y*');
grid on
title('B-PAM with E_b/N_0 = 5dB')
awgn_bpam_mod_15dB = awgn(bpam_mod_3e, 15+3);
scatterplot((A/2)*awgn_bpam_mod_15dB, 1, 0, 'y*');
grid on
title('B-PAM with E_b/N_0 = 15dB')


% 3.st
bit_size = 2;
EbN0 = 0:1:15;
[BER] = berawgn(EbN0, 'pam', bit_size);

n_bits_3st = 1000000;
k_3st = log2(bit_size);
snr_dB = EbN0+3+log10(k_3st);
y_noisy = zeros(n_bits_3st,length(snr_dB));
z_snr = zeros(n_bits_3st,length(snr_dB));
errVec = zeros(3,length(EbN0));
errcalc = comm.ErrorRate;
random_bits_3st = randi([0,1],[1000000,1]);
y_random_bits_pam = pammod(random_bits_3st, bit_size);

for jj = 1:length(snr_dB)
    reset(errcalc)
    y_noisy(:,jj) = awgn(real(y_random_bits_pam), snr_dB(jj),'measured');
    z_snr(:,jj) = pamdemod(complex(y_noisy(:,jj)),bit_size);
    errVec(:,jj) = errcalc(random_bits_3st, z_snr(:,jj));
end

figure (29)
semilogy(EbN0, BER, 'r');
hold on
semilogy(EbN0,errVec(1,:),'b.');
grid on
legend('Theoretical BER', 'Experimental BER')
title('Theoretical-Experimental BER for B-PAM')
xlabel('E_b/N_0')
ylabel('Bit Error Rate')
hold off


%--------------------------------------------------------------------------

% 4o Erwthma

% 4.a
qpsk_4 = comm.QPSKModulator(BitInput=true);
qpsk_bits_4 = A*qpsk_4(random_bits);
scatterplot(0.5*qpsk_bits_4, 1, 0, 'y*');
grid on
title('QPSK Modulated Signal')
text(0.5*A/sqrt(2)-0.12, 0.5*A/sqrt(2)-0.22, '00', 'Color', [1 1 1]);
text(-0.5*A/sqrt(2)-0.12, 0.5*A/sqrt(2)-0.22, '01', 'Color', [1 1 1]);
text(-0.5*A/sqrt(2)-0.12, -0.5*A/sqrt(2)+0.22, '11', 'Color', [1 1 1]);
text(0.5*A/sqrt(2)-0.12, -0.5*A/sqrt(2)+0.22, '10', 'Color', [1 1 1]);


% 4.b
awgn_4b_5db = awgn(qpsk_bits_4, 5);
scatterplot(awgn_4b_5db, 1, 0, 'y*');
grid on
title('QPSK with E_b/N_0 = 5dB')

awgn_4b_15db = awgn(qpsk_bits_4, 15);
scatterplot(awgn_4b_15db, 1, 0, 'y*');
grid on
title('QPSK with E_b/N_0 = 15dB')


% 4.c
bit_size_qpsk = 4;
EbN0 = 0:1:15;
[BER_4c] = berawgn(EbN0,'psk',bit_size_qpsk ,'diff');
figure(33)
semilogy(EbN0,BER_4c,'b')
grid on
title('Comparing BER for QPSK and BPSK modulation')
xlabel('E_b/N_0')
ylabel('Bit Error Rate')

bit_size_bpsk = 2;
n_bits_bpsk = 1000000;
k_bits_bpsk = log2(bit_size_bpsk);
snr_random_bits_bpsk = EbN0+10*log10(k_bits_bpsk);
y_noisy_bpsk = zeros(n_bits_bpsk,length(snr_random_bits_bpsk));
z_snr_bpsk = zeros(n_bits_bpsk,length(snr_random_bits_bpsk));
errVec_bpsk = zeros(3,length(EbN0));

errcalc = comm.ErrorRate;

x_random_bits_bpsk = randi([0 1],1000000,1);
y_random_bits_bpsk = pskmod(x_random_bits_bpsk, bit_size_bpsk);
for jj = 1:length(snr_random_bits_bpsk)
    reset(errcalc)
    y_noisy_bpsk(:,jj) = awgn(y_random_bits_bpsk,snr_random_bits_bpsk(jj),'measured');
    z_snr_bpsk(:,jj) = pskdemod(complex(y_noisy_bpsk(:,jj)), bit_size_bpsk);
    errVec_bpsk(:,jj) = errcalc(x_random_bits_bpsk,z_snr_bpsk(:,jj));
end

hold on
semilogy(EbN0,errVec_bpsk(1,:),'g.');

n_bits_qpsk = 1000000;
k_bits_qpsk = log2(bit_size_qpsk);
snr_random_bits_qpsk = EbN0+10*log10(k_bits_qpsk);
y_noisy_qpsk = zeros(n_bits_qpsk,length(snr_random_bits_qpsk));
z_snr_qpsk = zeros(n_bits_qpsk,length(snr_random_bits_qpsk));
errVec_qpsk = zeros(3,length(EbN0));

errcalc = comm.ErrorRate;

x_random_bits_qpsk = randi([0 1],1000000,1);
y_random_bits_qpsk = pskmod(x_random_bits_qpsk, bit_size_qpsk);

for jj = 1:length(snr_random_bits_qpsk)
    reset(errcalc)
    y_noisy_qpsk(:,jj) = awgn(y_random_bits_qpsk,snr_random_bits_qpsk(jj),'measured');
    z_snr_qpsk(:,jj) = pskdemod(complex(y_noisy_qpsk(:,jj)), bit_size_qpsk);
    errVec_qpsk(:,jj) = errcalc(x_random_bits_qpsk,z_snr_qpsk(:,jj));
end

semilogy(EbN0,errVec_qpsk(1,:),'r.');
legend('Theoretical BER','BPSK Empirical BER','QPSK Empirical BER');
hold off


% 4.d.i
text_file = fopen('rice_even.txt', 'r');
text_content = fscanf(text_file, '%c');
text_bits = dec2bin(text_content);

% 4.d.ii
bit_txt_dec = bin2dec(text_bits);
N_txt = 8;
L_txt = 2^N_txt ;
D_txt = (max(bit_txt_dec) - min(bit_txt_dec))/L_txt;
I_txt = round((bit_txt_dec - min(bit_txt_dec))/D_txt);
txt_quantized = min(bit_txt_dec) + I_txt*D_txt;
figure(35);
n_txt = 0:1:498;
stem (n_txt,txt_quantized,'rx-');
xlabel('Characters')
ylabel('Decimal Representation of Ascii Code')
title('Quantized Text Signal')


% 4.d.iii
txt_quantized_bits = dec2bin(round(txt_quantized)) - '0';
quantized_bit_stream = reshape(txt_quantized_bits.', [], 1);
quantized_bit_stream(end+1) = 0;

txt_qpsk_mod = qpsk_4(quantized_bit_stream);
scatterplot(txt_qpsk_mod, 1, 0, 'y*');
grid on
title('Text QPSK Constellation Diagram')

% 4.d.iv
txt_qpsk_mod_noisy5dB = awgn(txt_qpsk_mod, 5);
txt_qpsk_mod_noisy15dB = awgn(txt_qpsk_mod, 15);

% 4.d.v
qpskdemod_txt = comm.QPSKDemodulator(BitOutput=true);
txt_qpsk_demod_5dB_dec = qpskdemod_txt(txt_qpsk_mod_noisy5dB);
txt_qpsk_demod_15dB_dec = qpskdemod_txt(txt_qpsk_mod_noisy15dB);
scatterplot(txt_qpsk_mod_noisy5dB);
scatterplot(txt_qpsk_mod_noisy15dB);

% 4.d.vi
[num_qpsk_demod_5dB, error_qpsk_demod_5dB] = symerr(quantized_bit_stream, txt_qpsk_demod_5dB_dec);
[num_qpsk_demod_15dB, error_qpsk_demod_15dB] = symerr(quantized_bit_stream, txt_qpsk_demod_15dB_dec);
ber_theoretical_5dB = qfunc(sqrt(10^(5/10)));
% ber_theoretical_5dB = 0.0377
ber_theoretical_15dB = qfunc(sqrt(10^(15/10)));
% ber_theoretical_15dB = 0

% 4.d.vii
txt_demod_5dB  = char(bin2dec(reshape(char('0' + txt_qpsk_demod_5dB_dec(1:3395)),7,[]).'))';
txt_demod_15dB  = char(bin2dec(reshape(char('0' + txt_qpsk_demod_15dB_dec(1:3395)),7,[]).'))';
txt_5dB = fopen('txt_5db_03119908.txt', 'w');
fprintf(txt_5dB, '%c', txt_demod_5dB);
txt_15dB = fopen('txt_15db_03119908.txt', 'w');
fprintf(txt_15dB, '%c', txt_demod_15dB);


%------------------------------------------------------------

% 5o Erwthma

% 5.a
[sound_wave, fs_sound] = audioread('soundfile2_lab2.wav'); 
dt_sound = 1/fs_sound;
t_sound = 0:dt_sound:(length(sound_wave)*dt_sound)-dt_sound;
figure (39);
plot(t_sound, sound_wave);
grid on
xlabel('Time [s]')
ylabel('Amplitude')
title('Sound Signal')


% 5.b
sound_bits = 8;
sound_levels = 2^sound_bits;
sound_step_size = (max(sound_wave) - min(sound_wave))/sound_levels;
sound_index = round((sound_wave - min(sound_wave))/sound_step_size);
quantized_sound = min(sound_wave) + sound_index * sound_step_size;
figure (40);
stem(t_sound,quantized_sound, 'gx-');
grid on
xlabel('Time [s]')
ylabel('Amplitude')
title('Quantized Sound Signal at 8 bits')


% 5.c
quantized_sound_mult = quantized_sound * 1000;
quantized_sound_bits = dec2bin(quantized_sound_mult,16) - '0';
quantized_sound_bit_stream = reshape(quantized_sound_bits.',[],1);
qpsk_mod_5 = comm.QPSKModulator('BitInput', true);
sound_qpsk_mod = qpsk_mod_5(quantized_sound_bit_stream);
scatterplot(sound_qpsk_mod, 1, 0, 'y*');
grid on
title('Sound QPSK Constellation Diagram')


% 5.d
awgn_sound_qpsk_mod_4dB = awgn(sound_qpsk_mod, 4);
awgn_sound_qpsk_mod_14dB = awgn(sound_qpsk_mod, 14);


% 5.e
qpsk_demod_5 = comm.QPSKDemodulator('BitOutput',true); 
sound_qpsk_demod_4dB = qpsk_demod_5(awgn_sound_qpsk_mod_4dB);
sound_qpsk_demod_14dB = qpsk_demod_5(awgn_sound_qpsk_mod_14dB);
scatterplot(awgn_sound_qpsk_mod_4dB);
grid on
title('Sound QPSK with E_s/N_0 = 4dB')
scatterplot(awgn_sound_qpsk_mod_14dB);
grid on
title('Sound QPSK with E_s/N_0 = 14dB')


% 5.st
[num_qpsk_demod_4dB, error_qpsk_demod_4db] = symerr(quantized_sound_bit_stream, sound_qpsk_demod_4dB);
[num_qpsk_demod_14dB, error_qpsk_demod_14db] = symerr(quantized_sound_bit_stream, sound_qpsk_demod_14dB);
ber_theoretical_4dB = qfunc(sqrt(10^(4/10)));
% ber_theoretical_4dB = 0.0565
ber_theoretical_14dB = qfunc(sqrt(10^(14/10)));
% ber_theoretical_14dB = 0


% 5.z
sound_reshaped_14dB = char('0'+ (reshape(sound_qpsk_demod_14dB, 16, []).'));
sound_demod_dec_14dB = typecast(uint16(bin2dec(sound_reshaped_14dB)),'int16');
sound_demod_14dB = double(abs(sound_demod_dec_14dB))./1000 ;

sound_reshaped_4dB = char('0'+ (reshape(sound_qpsk_demod_4dB, 16, []).'));
sound_demod_dec_4dB = typecast(uint16(bin2dec(sound_reshaped_4dB)),'int16');
sound_demod_4dB = double(abs(sound_demod_dec_4dB))./1000 ;

sound_demod_14dB_norm = sound_demod_14dB / max(sound_demod_14dB);
audiowrite("sound_14dB_03119908.wav", sound_demod_14dB_norm , 44100)
sound_demod_4dB_norm = sound_demod_4dB / max(sound_demod_4dB);
audiowrite("sound_4dB_03119908.wav", sound_demod_4dB_norm, 44100)

%}




