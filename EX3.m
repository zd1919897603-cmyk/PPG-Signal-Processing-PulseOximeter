clc
close all
clear 

%%
% Step 1
data = readmatrix('pulse.txt');
Led_R  = data(:,1);
Led_IR = data(:,2);

%%
% Step 2

fs = 100;           
T_offset = 10;      
T_duration = 60;    

start_index = T_offset * fs + 1;

end_index = start_index + (T_duration * fs) - 1; 


if end_index > length(Led_R)
    end_index = length(Led_R);
    warning('Data is shorter than 70 seconds, capturing until the end.');
end

Led_R_valid = Led_R(start_index : end_index);
Led_IR_valid = Led_IR(start_index : end_index);

t = (0:length(Led_R_valid)-1)/fs;

figure('Color','w');hold on;
subplot(2,1,1);
plot(t,Led_R_valid,'LineWidth',0.5);
title('RED');
xlabel('Time[s]');
xlim([0,60]);ylim([2.004e5,2.016e5]);

subplot(2,1,2);
plot(t,Led_IR_valid,'LineWidth',0.5);
title('INFRARED');
xlabel('Time[s]');
xlim([0,60]);ylim([2.62e5,2.65e5]);

%%
% Step 3 low_pass filter
clear pi
fp = 3;
fs_stop = 6;
Rp = 1;
As = 60;
fs_sample = 100;

Wp = 2*pi*fp;
Ws = 2*pi*fs_stop;

N_min = log10((10^(As/10)-1)/(10^(Rp/10)-1)) / (2 * log10(Ws/Wp));
N = ceil(N_min);

Wc1 = Wp/((10^(Rp/10)-1)^(1/(2*N)));
Wc2 = Ws/((10^(As/10)-1)^(1/(2*N)));

Wc_avg = (Wc1 + Wc2)/2;
fc = round(Wc_avg/(2*pi),3);
Wc = 2*pi*fc;

k = 0:N-1;
poles = Wc*exp(1j*(pi/2+(2*k+1)*pi/(2*N)));

%Ha[s] 
[num_s, den_s] = zp2tf([],poles',Wc^N);
figure('Name','left half-plane s');
zplane([],poles');title('Poles(left half-plane)');

f_plot = 0:0.1:40;
W_plot = 2*pi*f_plot;
Ha_mag2 = 1./(1 + (W_plot./ Wc).^(2*N));
figure; plot(f_plot, Ha_mag2); grid on;
title('The squared magnitude of the ananlog filter'); xlabel('Frequency (Hz)'); ylabel('|Ha(j\Omega)|^2');

[num_z, den_z] = impinvar(num_s, den_s,fs);
figure;
freqz(num_z, den_z,1024,fs);
title('The digital filter');

Led_R_filtered = filtfilt(num_z, den_z,Led_R_valid);
Led_IR_filtered = filtfilt(num_z, den_z,Led_IR_valid);

% figure('Color','w');
% subplot(2,1,1);
% plot(Led_R_valid,'Color',[0 1 0]);hold on;
% plot(Led_R_filtered,'r','LineWidth',1);
% title('Red-Original vs Filterd');legend('Original','Filtered');
% 
% subplot(2,1,2);
% plot(Led_IR_valid,'Color',[1 0.6 0]);hold on;
% plot(Led_IR_filtered,'b','LineWidth',1);
% title('Ifrared-Original vs Filterd');legend('Original','Filtered');

%%
% Step 4 High_pass filter
fs = 100;
fst = 0.05;
fp = 0.75;
Wp = 2*pi*fp;
Ws = 2*pi*fst;
fc = (fst + fp)/(2*fs);
df = (fp - fst)/fs;
N = ceil(5.5/df);
M = ceil((N -1)/2);
n = 0:1:N-1;

h_id = 2*fc.*sinc(2*fc.*(n - M));
h_id = -h_id;
h_id(M+1) = h_id(M+1)+1;
Wn = 0.42 - 0.5.*cos((2*pi.*n)./(N -1)) + 0.08.*cos((4*pi.*n)./(N -1));
hn = h_id .* Wn;

[h, f] = freqz(hn,1, 1024,fs);
figure;
plot(f,20*log10(abs(h)),'LineWidth',1);
grid on;hold on;
yline(-60, '--r', 'Stop-band -60dB');
xline(fst, '--g', 'fst=0.05');
xline(fp, '--b', 'fp=0.75');
title('Step 4: FIR High-pass Filter Design');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
ylim([-100, 10]);

Led_R_final = filtfilt(hn,1,Led_R_filtered);
Led_IR_final = filtfilt(hn,1,Led_IR_filtered);

%%
%Step 5- plot the filtered signals
figure;
subplot(3,1,1);
plot(t,Led_R_valid,'b','LineWidth',0.5);
title('original Signal');
xlim([0,60]);ylim([2.005e5,2.015e5]);
xlabel('Time[s]');
subplot(3,1,2);
plot(t,Led_R_filtered,'b','LineWidth',0.5);
title('after low pass');
xlim([0,60]);ylim([2.005e5,2.015e5]);
xlabel('Time[s]');
subplot(3,1,3);
plot(t,Led_R_final,'b','LineWidth',0.5);
title('after high pass');
xlim([0,60]);ylim([-200,200]);
xlabel('Time[s]');

figure
subplot(3,1,1);
plot(t,Led_IR_valid,'b','LineWidth',0.5);
title('original InfraRed Signal');
xlim([0,60]);ylim([2.62e5,2.65e5]);
xlabel('Time[s]');
subplot(3,1,2);
plot(t,Led_IR_filtered,'b','LineWidth',0.5);
title('after low pass');
xlim([0,60]);ylim([2.62e5,2.65e5]);
xlabel('Time[s]');
subplot(3,1,3);
plot(t,Led_IR_final,'b','LineWidth',0.5);
title('after high pass');
xlim([0,60]);ylim([-500,500]);
xlabel('Time[s]');

%%
%Step 6 - Pulse rate computation
L = length(Led_R_final);
N = 2^nextpow2(L);

Y = fftshift(fft(Led_R_final,N));
Y_dB = 20*log10(abs(Y));
f2 = (-N/2 : N/2-1) * (fs/N);
[max_pulse_rate,max_idx] = max(Y_dB);
f_pulse_rate = f2(max_idx);
f_pulse_rate = -f_pulse_rate;

figure;
plot(f2,Y_dB,'LineWidth',1);
title('Filtered red signal-BPM =',60*f_pulse_rate);
xlim([-50,50]);ylim([20,110]);
xlabel('Frequency[Hz]');ylabel('Amplitude[dB]');
xline(max_pulse_rate,'r','BPM');

%%
%Step 7 - Saturation computation
% Red signal
[pks_R_ac, locs_p_R_ac] = findpeaks(Led_R_final, 'MinPeakDistance', 60);
[vls_R_ac, locs_v_R_ac] = findpeaks(-Led_R_final, 'MinPeakDistance', 60);
vls_R_ac = -vls_R_ac; 

[vls_R_dc, locs_v_R_dc] = findpeaks(-Led_R_filtered, 'MinPeakDistance', 60);
vls_R_dc = -vls_R_dc;
t3 = (0:length(Led_R_final)-1)/fs;
Upper_R_ac = interp1(locs_p_R_ac/fs, pks_R_ac, t3, 'spline');
Lower_R_ac = interp1(locs_v_R_ac/fs, vls_R_ac, t3, 'spline');
IDC_R_curve = interp1(locs_v_R_dc/fs, vls_R_dc, t3, 'spline');

IAC_R_curve = Upper_R_ac - Lower_R_ac; 

% IR signal
% AC part: Led_IR_final
[pks_IR_ac, locs_p_IR_ac] = findpeaks(Led_IR_final, 'MinPeakDistance', 60);
[vls_IR_ac, locs_v_IR_ac] = findpeaks(-Led_IR_final, 'MinPeakDistance', 60);
vls_IR_ac = -vls_IR_ac; 

% DC part: Led_IR_lp 
[vls_IR_dc, locs_v_IR_dc] = findpeaks(-Led_IR_filtered, 'MinPeakDistance', 60);
vls_IR_dc = -vls_IR_dc;

t4 = (0:length(Led_IR_final)-1)/fs;

% Spline Interpolation
Upper_IR_ac = interp1(locs_p_IR_ac/fs, pks_IR_ac, t4, 'spline');
Lower_IR_ac = interp1(locs_v_IR_ac/fs, vls_IR_ac, t4, 'spline');
IDC_IR_curve = interp1(locs_v_IR_dc/fs, vls_IR_dc, t4, 'spline');

% Peak to peak
IAC_IR_curve = Upper_IR_ac - Lower_IR_ac;

R_signal = (IAC_R_curve ./ IDC_R_curve) ./ (IAC_IR_curve ./ IDC_IR_curve);
R_mean = mean(R_signal);
SaO2 = 110 - 25 * R_mean;

%%(Final Plots)

figure;

subplot(2,2,1);
plot(t3, Led_R_final, 'k','LineWidth',0.1); hold on;
plot(t3, Upper_R_ac, '-ro', 'MarkerSize', 1);
plot(t3, Lower_R_ac, '-ro','MarkerSize',1);   
grid on;
title('RED Signal: AC Envelopes');
xlabel('Time [s]'); ylabel('Amplitude');

subplot(2,2,2);
plot(t3, Led_R_filtered, 'k','LineWidth',0.1); hold on;
plot(t3, IDC_R_curve, 'r', 'LineWidth', 1.2);
grid on;
title('RED Signal: DC Component (IDC)');
xlabel('Time [s]'); ylabel('Intensity');

subplot(2,2,3);
plot(t4, Led_IR_final, 'k','LineWidth',0.1); hold on;
plot(t4, Upper_IR_ac, 'r', 'LineWidth', 1);
plot(t4, Lower_IR_ac, 'r', 'LineWidth', 1);
grid on;
title('IR Signal: AC Envelopes');
xlabel('Time [s]'); ylabel('Amplitude');


subplot(2,2,4);
plot(t4, Led_IR_filtered, 'k', 'LineWidth', 0.1);hold on;
plot(t4,IDC_IR_curve,'r','LineWidth',1)
grid on;
title('IFRARED Signal: DC Component (IDC)');
xlabel('Time [s]'); ylabel('Intensity');

sgtitle(['SO2 =',num2str(SaO2, '%.2f'), '%']);