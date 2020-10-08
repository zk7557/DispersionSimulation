%% initial setup
% um, fs unit
global c;
c = 0.3;

% time
tspan = 1000;
dt = 0.1;
t = -tspan:dt:tspan;
ft = fftx(t);

%% create a pulse
pulseWidth = 300; % fs
centralWavelength = 1; % um
centralFrequency = c / centralWavelength;

% gaussian pulse
% pulse = 1/pulseWidth / sqrt(2*pi)* exp(- t.^2/(2*pulseWidth^2)) .* ...
% 		  exp(1i * 2*pi*c/centralWavelength * t);

% triple harmonic pulse, square window
% t1 = -tspan:dt:-tspan*0.1;
% t2 = -tspan*0.1+dt:dt:tspan*0.1-dt;
% t3 = tspan*0.1:dt:tspan;
% [~, nt] = size(t1);
% pulse = sin(2*pi*0.3*t2)+sin(2*pi*0.3*t2*3)+sin(2*pi*0.3*t2*5);
% pulse = [zeros(1,nt) pulse zeros(1,nt)];
% clearvars t1 t2 t3 nt;

% double harmonic complex pulse, "Fermi-Dirac" window
buffer = 20;
amplitude = 1./(1 + exp((abs(t)-pulseWidth) ./ (buffer)));
pulse = amplitude .* (exp(2i*pi* centralFrequency*t) + ...
                      exp(2i*pi* centralFrequency*t * 2));
clearvars buffer amplitude;



% check the shape of the pulse
figure(1); plot(t, real(pulse)); title('initial pulse diagram');

% a nice plot of spectrum as function of time
instSptr = instantSpectrum(pulse, t, 5);
figure(2);
spectrum = imagesc(instSptr, ...
		   'xData', -tspan:5:tspan, 'yData', ft);
colorbar; title('Initial Spectrum');

% peak frequency vs time
%figure(5);


%% apply temporal dispersion to the pulse
fpulse = fft(pulse);
shiftedPulse = shiftedIFFT(fpulse, ft);

% plot to check the pulse shape
figure(3); plot(t, real(shiftedPulse));
title('Pulse after Dispersion');
xlabel('time (fs)'); ylabel('real E (au)');

% also instant spectrum
figure(4);
spectrum2 = imagesc(instantSpectrum(shiftedPulse, t, 5), ...
		   'xData', -tspan:5:tspan, 'yData', ft);
colorbar; title('Spectrum after Disperison');


%% find the GDD and TOD of a configuration
global c;
c = 0.3;
df = 0.0005;
f = 0:df:1-df;  
cf = 0.3;        % central frequency

% calculate group delay
% cache: gratingDelay(1, 1, pi/6, f, cf)
% cache: prismDelay(pi/12, 10000, 0, f, cf)
tg = prismDelay(pi/3, 20000, pi/6, f, cf);

% plot group delay
figure(21); plot(f,tg); axis([0 1 -2*10^4 10^5]);


% calculate GDD and TOD
GDD = (tg(2:end) - tg(1:end-1))/(df*2*pi);
TOD = (GDD(2:end) - GDD(1:end-1))/(df*2*pi);

% show the GDD/TOD ratio
ratio23 = GDD(600)/TOD(600);

% plot GDD and TOD
rng = (550:650);
figure(22); plot(f(rng),GDD(rng));
figure(23); plot(f(rng),TOD(rng));

%% Test area
global c;
prismDelay(1, 2, 3, 4, 4)