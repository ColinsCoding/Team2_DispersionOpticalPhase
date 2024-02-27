clear all;
close all;
clc;




D1 = -600;
D2 = -900;
distance = 5/1000;



sample_rate = 200e9;
shift_freq = 40e9;
time_step = 1/sample_rate;
%absorb_spread = 0.5e-9;
absorb_spread = 0.5e-10;


time_length = 50e-9;

time = linspace(-10e-9,time_length,fix(2*time_length*sample_rate)+1);
freq = linspace(0, sample_rate, length(time));


T0 = 20e-10;
time_delay_of_pulse = 20e-9;
chirp_mag = exp(-0.5.*(1).*((time-time_delay_of_pulse)/T0).^2);
for i = 1:(length(time)-1)/2
    if time(i) <= 0
        chirp_mag(i) = 0;
    end
end
chirp_phase = exp(-0.5.*(0.3i).*((time-time_delay_of_pulse)/T0).^3);

chirp_mag = chirp_mag.*chirp_phase;

Dw1 = exp(1i*(D1/(2*pi*299792)).*distance.*freq.*freq);
Dw2 = exp(1i*(D2/(2*pi*299792)).*distance.*freq.*freq);
Dw1_inv = log(1i*(D1/(2*pi*299792)).*freq.*freq);
Dw2_inv = log(1i*(D2/(2*pi*299792)).*freq.*freq);



chirp_ft = fft(chirp_mag);

after_D1 = chirp_ft.*Dw1;
after_D2 = chirp_ft.*Dw2;

D1_time = ifft(after_D1);
D2_time = ifft(after_D2);

D1_time_intensity = D1_time.*conj(D1_time);
D2_time_intensity = D2_time.*conj(D2_time);




for i = 1:length(time)
    phase(i) = -180 + 360*rand();
end

phase1 = phase.*pi/180;
phase2 = phase1;


error = [];
itt_cnt = 0;
itt_num = 1000;

f1_meas = sqrt(D1_time_intensity);
f2_meas = sqrt(D2_time_intensity);


% f = figure;
% p = uipanel(f,"Position",[0.1 0.1 0.8 0.8],...
%     "BackgroundColor","w");
% ax = axes(p);
% 
% myVideo = VideoWriter('myVideoFile'); %open video file
% myVideo.FrameRate = 30;  %can adjust this, 5 - 10 works well for me
% open(myVideo)
% tic
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

%     plot(time/1e-9, ifft(fft(f1_meas.*exp(1i.*phase1)).*Dw2./Dw1).*conj(ifft(fft(f1_meas.*exp(1i.*phase1)).*Dw2./Dw1)),...
%         time/1e-9, smoothdata(ifft(fft(f1_meas.*exp(1i.*phase1)).*Dw2./Dw1).*conj(ifft(fft(f1_meas.*exp(1i.*phase1)).*Dw2./Dw1)),'sgolay'),...
%         time/1e-9,f2_meas.*conj(f2_meas))
%     legend(["GSA Phase" "Smoothed GSA Phase" "Target"])
%     u.Value = i;
%     M(i) = getframe(gcf);
%     writeVideo(myVideo, M(i));
%     ylim([0 0.5])
%     %pause(0.1);

end
%toc
% close all
% close(myVideo)
plot(freq, 10*log10(abs(fft(f1_meas.*exp(1i.*phase1)))))




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
% plot(time/1e-9, f2_meas.*conj(f2_meas),"-","LineWidth",2, Color=[0 0 0.6]  )
% hold on
% plot(time/1e-9, smoothdata(ifft(fft(f1_meas.*exp(1i.*phase1)).*Dw2./Dw1).*conj(ifft(fft(f1_meas.*exp(1i.*phase1)).*Dw2./Dw1)),'sgolay'), "--", "LineWidth",2, Color=[1 0 1])
% plot(time/1e-9, ifft(fft(f1_meas).*Dw2./Dw1).*conj(ifft(fft(f1_meas).*Dw2./Dw1)), "r-.", LineWidth=2)
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
% subplot(1,2,1)
% plot(time/1e-9, abs(chirp_mag), LineWidth=2)
% xlabel("time (ns)")
% ylabel("Intensity (a.u.)")
% title("Chirp Intensity")
% 
% subplot(1,2,2)
% plot(time/1e-9, unwrap(angle(chirp_mag)), LineWidth=2)
% xlabel("time (ns)")
% ylabel("Phase (radians)")
% title("Chirp Phase")



