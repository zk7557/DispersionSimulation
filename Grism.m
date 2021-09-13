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
% Our gratings: gratingDelay(1, 10000, pi/6, f, cf)
% Equilateral BK7: prismDelay(pi/3, 20000, 0.507, f, cf)
% Thorlabs AFS-FS: prismDelay(69.1* pi/180, 20000,pi/4, f, cf)
tg = gratingDelay(1, 10000, pi/6, f, cf);

% plot group delay
figure(21); plot(f,tg); axis([0 1 -2*10^4 10^5]);


% calculate GDD and TOD
GDD = (tg(2:end) - tg(1:end-1))/(df*2*pi);
TOD = (GDD(2:end) - GDD(1:end-1))/(df*2*pi);

% show the GDD/TOD ratio
ratio23 = GDD(600)/TOD(600);

% plot GDD and TOD
rng = (550:650);
figure(22); plot(f(rng),GDD(rng));title('GDD');
figure(23); plot(f(rng),TOD(rng));title('TOD');

%% Test area
global c;
prismDelay(1, 2, 3, 4, 4)

%% functions
% refractive index for BK7, function of lambda
function n = rIndexL(wavelength)
	w2 = wavelength .^ 2;
	n = sqrt(1 + 1.03961212*w2./(w2-0.00600069867) + ...
				 0.231792344*w2./(w2-0.0200179144) + ...
				 1.01046945*w2./(w2-103.560653) ...
			);
end

% refractive index for BK7, function of frequency
function n = rIndexF(f)
	global c;
	wavelength = c ./ f;
	n = rIndexL(wavelength);
end

% instant spectrum, gaussian window
function s = instantSpectrum(wave, t, window)
	dt = t(2) - t(1);
	[~, nRow] = size(t);
	[~, sizeWave] = size(wave);
	if nRow ~= sizeWave
		error('wave - t dimention mismatch');
	end
	nCol = floor(nRow / (window / dt));
	
	% filtered wave
	fWave = zeros(nRow, nCol);
	for i = 1:nCol
		filter = 1/sqrt(2*pi)/window * ...
				 exp(-( t - ((i-1)*window+t(1)) ).^2 / 2 / window^2);
		fWave(:,i) = wave .* filter;
	end
	s =abs(fft(fWave));
end

% frequency axes of fft
function f = fftx(t)
	[~, n] = size(t);
	dt = t(2) - t(1);
	df = 1/(dt * n);
	f = 0:df:df*(n-1);
end

% time delay from grism in Louradour_2012_OptExp
% arguments: apex angle, gratings pitch, grism distance T2O, incident angle
%			 light frequency as array, the zero-dispersion primary freq
function tg = grismDelay(a, d, dd, theta, f)
	global c;
	n = rIndexF(f);
	%lambda over d
	ld = c/d./f;
	% angle out of the first grism
	theta5 = asin(n .* sin(pi/2 - a - asin( ...
			 sin(a - asin(sin(theta) / n)) - ld )));
	% path length between the grisms
	c1c2 = dd ./ cos(theta5);
	
	% the common denomitator of cba and ai2
	usefulAngle = a + asec(csc(theta)*n);
	denominator =	sqrt(1- (2*pi*c + 2*pi*f*d *cos( ...
					usefulAngle)).^2./(2*pi*d*f).^2) .* ...
					sqrt(1 - (sin(theta)/n).^2);
	% C2B2A2
	xp = 1; % this should be cancelled when taking the differential
	cba =	(xp * sin(a) * (sqrt(1 - (sin(theta)/n).^2) + ...
			sin(a - asin(ld + cos(usefulAngle))))) ...
			./ denominator;
	% A2I2
	ai2 = - (xp * cos(a + acos(ld + cos(usefulAngle))) ...
			* sin(usefulAngle)) ./ denominator;
	
	% total time delay
	tg = 2/c * (n .* cba + c1c2 - ai2);
	
end

% time delay from prism compressor, shown in CompressorAnnotation.pdf
% a:    apex angle
function tp = prismDelay(a, dd, theta, f, pf)
    global c;
	f = [f pf];
	n = rIndexF(f);
	n2 = n .^2;
	
	xp = 1; % useless variable - well, not that useless

    theta3 = asin(n .* sin(a - asin(sin(theta) ./ n)));
    cos1 = cos(asin(sin(theta) ./ n));
    
    tp = 1/c * ( ...
                dd ./ cos(theta3) + ...
                n .* (xp - dd * sin(a) * tan(theta3) ./ cos1) - ...
                xp + dd * tan(theta3) .* ...
                cos(a - asin(sin(theta) ./ n))./ cos1 .* sin(theta) ...
                );
    
    tp(1) = 0;
    tp = real(tp - tp(end));
    tp = tp(1:end-1);
end

% time delay from grating compressor, the one in use
% d:    groove distance, in um
% dd:   distance between the gratings, in um
% theta:incident angle, respact to normal direction
% f:    frequency array
% pf:   central frequency, only used as zero reference
function tc = gratingDelay(d, dd, theta, f, pf)
	global c;
    f = [f pf];
	sinThetaPrime = sin(theta) - (c / d) ./  f;
	tc = cos(theta - asin(sinThetaPrime)) + ...
		 dd ./ sqrt(1 - sinThetaPrime .^2);
    tc = real(tc - tc(end));
    tc(1)=0;
    tc = tc(1:end-1);
end

% custom fourier transform
function X = myFFT(x)
	[~, N] = size(x);
	% exp(2 pi i * t * f), t is the integrand, f is the conjugate
	M = exp(2i*pi/N * (0:(N-1)) .* (0:(N-1))' );
	X = M * x';
	X = X';     % output as row array
	% identical to fft, now we can modify it
end

function X = shiftedIFFT(x, ft)
	[~, N] = size(x);
	% exp(-2 pi i * f * t), f is the integrand, t is the conjugate
	
	% three matrices, f, t, delay(f), horizontal f first
	Mf = (0:N-1) .* ones(N,1);
	Mt = ones(1,N) .* (0:N-1)';
    
    % delay function of frequency, use the following format:
    % Mdelay = delayFunction(ft) .* ones(N, 1); 
    % where delayfunction returns a 1xN row vector
	Mdelay = prismDelay(69.1* pi/180, 40000 ,pi/4, ft, 0.3) .* ones(N,1) ;
    
	M = Mf .*(Mt - Mdelay);
	M = exp((-2i*pi/N) * M)/N;
	
	X = M * x';
	X = X';     % output as row array
	
end



