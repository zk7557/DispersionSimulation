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
    
%     tp(1) = 0;
    tp = real(tp - tp(end));
    tp = tp(1:end-1);
end

% time delay from grating compressor, the one in use
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
	Mdelay = prismDelay(pi/12, 10000, 0.1, ft, 0.3) .* ones(N,1) ;
    
	M = Mf .*(Mt - Mdelay);
	M = exp((-2i*pi/N) * M)/N;
	
	X = M * x';
	X = X';     % output as row array
	
end