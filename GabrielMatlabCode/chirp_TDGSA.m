clear all; 
%close all; 
clc;


D1 = -353;
D2 = -872;
distance = 5/1000;


sample_rate = 200e9;
shift_freq = 40e9;
time_step = 1/sample_rate;
absorb_spread = 0.2e-9;


time_length = 50e-9;

time = linspace(-10e-9,time_length,fix(2*time_length*sample_rate)+1);
freq = linspace(0, sample_rate, length(time));


C = 0;
T0 = 20e-10;
time_delay_of_pulse = 20e-9;
chirp = exp(-0.5.*(1+1i*C).*((time-time_delay_of_pulse)/T0).^2);
% for i = 1:(length(time)-1)/2
%     if time(i) <= 0
%         chirp(i) = 0;
%     end
% end


ampt_chirp = sqrt(abs(chirp));
phase_chirp = angle(chirp);

chirp_ft = fft(sqrt(chirp));

for i = 1:length(freq)
    absorb_func(i) = (absorb_spread*freq(i)-absorb_spread*shift_freq)^2/(1+(absorb_spread*freq(i)-absorb_spread*shift_freq)^2) + ...
        (absorb_spread*freq(i)-absorb_spread*shift_freq*2)^2/(1+(absorb_spread*freq(i)-absorb_spread*shift_freq*2)^2) + ...
        (absorb_spread*freq(i)-absorb_spread*shift_freq*0.5)^2/(1+(absorb_spread*freq(i)-absorb_spread*shift_freq*0.5)^2);
end

absorb_func = absorb_func - (min(absorb_func) - 10^-5);
absorb_func = absorb_func./max(absorb_func);

absorb_func = absorb_func.*exp(1i.*freq./absorb_func);


after_gas_chirp = chirp_ft.*absorb_func;

after_gas_time = ifft(after_gas_chirp);


Dw1 = exp(1i*(D1/(2*pi*299792)).*distance.*freq.^2);
Dw2 = exp(1i*(D2/(2*pi*299792)).*distance.*freq.^2);

% Dw1_inv = log(1i*(D1/(2*pi*299792)).*freq.^2);
% Dw2_inv = log(1i*(D2/(2*pi*299792)).*freq.^2);



after_disp_1 = after_gas_chirp.*Dw1;
after_disp_2 = after_gas_chirp.*Dw2;


ad1 = ifft(after_disp_1);
ad2 = ifft(after_disp_2);

for i = 1:length(time)
    phase(i) = -180 + 360*rand();
end

phase1 = phase.*pi/180;
error = [];
itt_num = 1000;


f1_int = ad1.*conj(ad1);
f2_int = ad2.*conj(ad2);

f1_meas = sqrt(f1_int);
f2_meas = sqrt(f2_int);
tic
for i = 1:itt_num
    f1_cmplx = f1_meas.*exp(1i.*phase1);
    f1_ft = fft(f1_cmplx);

    tx12 = f1_ft.*Dw2./Dw1;

    tx12_time = ifft(tx12);


    t2_phase = angle(tx12_time);

    f2_cmplx = f2_meas.*exp(1i.*t2_phase);

    f2_ft = fft(f2_cmplx);

    tx21 = f2_ft.*Dw1./Dw2;

    tx21_time = ifft(tx21);

    phase1 = angle(tx21_time);


    %error = [error; max(abs(f1_cmplx - tx21_time))];
    error = [error; max(abs(f1_cmplx.*conj(f1_cmplx) - tx21_time.*conj(tx21_time)))];
    error_percent(i) = max(abs((f1_cmplx.*conj(f1_cmplx) - tx21_time.*conj(tx21_time)).*100./(f1_cmplx.*conj(f1_cmplx))));


end

% for i = 1:itt_num
%     phase1 = angle(ifft(fft(f2_meas.*exp(1i.*angle(ifft(fft(f1_meas.*exp(1i.*phase1)).*Dw2./Dw1)))).*Dw1./Dw2));
% end

toc

f1_with_phase = f1_meas.*exp(1i.*phase1);

f1_fft = fft(f1_with_phase);
f1_fft_smooth = smoothdata(f1_fft, 'sgolay');

plot(freq, f1_fft_smooth)

%plot(time, (abs(ifft(after_gas_chirp.*Dw1.*exp(1i.*phase1)))), time, abs(f1_meas.*exp(1i.*phase1)))








% figure()
% plot(time/1e-9, f2_meas.*conj(f2_meas),"-","LineWidth",3, Color=[0 0 0.6]  )
% hold on
% plot(time/1e-9, ifft(fft(f1_meas.*exp(1i.*phase1)).*Dw2./Dw1).*conj(ifft(fft(f1_meas.*exp(1i.*phase1)).*Dw2./Dw1)), "--", "LineWidth",0.1, Color=[1 0 1])
% plot(time/1e-9, ifft(fft(f1_meas).*Dw2./Dw1).*conj(ifft(fft(f1_meas).*Dw2./Dw1)), 'r')
% legend(["Intensity measured from line 2" "Intensity from line 1 with TDGSA phase at line 2" "Intensity from line 1 without TDGSA phase at line 2"])
% xlabel("Time (ns)")
% ylabel("Tntensity (a.u.)")
% title("Intensity")
% 
% figure()
% plot(time/1e-9, smoothdata(f2_meas.*conj(f2_meas),'sgolay'),"-","LineWidth",2, Color=[0 0 0.6]  )
% hold on
% plot(time/1e-9, smoothdata(ifft(fft(f1_meas.*exp(1i.*phase1)).*Dw2./Dw1).*conj(ifft(fft(f1_meas.*exp(1i.*phase1)).*Dw2./Dw1)),'sgolay'), "--", "LineWidth",2, Color=[1 0 1])
% plot(time/1e-9, smoothdata(ifft(fft(f1_meas).*Dw2./Dw1).*conj(ifft(fft(f1_meas).*Dw2./Dw1)),'sgolay'), "r-.", LineWidth=2)
% legend(["Intensity measured from line 2" "Intensity from line 1 with TDGSA phase at line 2" "Intensity from line 1 without TDGSA phase at line 2"])
% xlabel("Time (ns)")
% ylabel("Tntensity (a.u.)")
% title("Intensity After Smoothing")
% 
% figure
% plot(error, LineWidth=2)
% xlabel("iteration number")
% ylabel('Max\{ABS(|F_{1}(t)|^{2} - |IFFT\{FFT\{F_{2}(t)\}D^{-1}_{2}(\omega)D_{1}(\omega)\}|^{2})\}')
% title("Error")
% 
% 
% figure()
% subplot(2,1,1)
% plot(freq/1e9, 10*log10(abs(absorb_func)), LineWidth=2)
% xlabel("Frequency (GHz)")
% ylabel("Absorbtion (dB)")
% title("Absorption spectrum Magnitude")
% 
% 
% subplot(2,1,2)
% plot(freq/1e9, unwrap(angle(absorb_func)), LineWidth=2)
% xlabel("Frequency (GHz)")
% ylabel("Phase (radian)")
% title("Absorption spectrum phase")
% 
% 
% 
% figure()
% subplot(1,2,1)
% plot(time/1e-9, abs(chirp), LineWidth=2)
% xlabel("time (ns)")
% ylabel("Intensity (a.u.)")
% title("Chirp Intensity")
% 
% subplot(1,2,2)
% plot(time/1e-9, unwrap(angle(chirp)), LineWidth=2)
% xlabel("time (ns)")
% ylabel("Phase (radians)")
% title("Chirp Phase")









