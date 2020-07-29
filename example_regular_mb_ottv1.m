% Generate a regular Mathieu beam using OTTv1 BscBessel basis
%
% This uses a method similar to section 1.3 from
%
%   Inverse-Fourier Non-diffracting Beams for Optical Trapping
%   Martinez-Ruiz et al. [Journal TBD], 2020.
%
% This functionality is planned for OTTv2/vswf.BesselBasis.MathieuBeam
% Copyright Isaac Lenton, 2020

% Nmax is the Nmax of the final beam.  The beam is only valid inside
% the corresponding Nmax sphere.  Lmax is related to the azimuthal
% resolution of the Mathieu beam ring.
Nmax = 50;        % VSWF mode number (truncation number)
Lmax = 100;      % Azimuthal mode number (truncation number)

morder = 11;       % m-order for beam
ellip = 40;      % Ellipticity of beam
parity = 'odd';  % Partity of beam (either even or odd)
theta = 0.5*pi/2;     % Angle for the non-diffracting ring in k-space (radians)

% For helical beams, create an even and odd beam and add them with a
% complex phase term.

%% Do point matching to determine modes

Npts = 2*Lmax+1;
phi = linspace(0, 2*pi, Npts+2);
phi = phi(1:end-1) + (phi(2) - phi(1))./2;

% Calculate MB
switch parity
  case 'odd'
    A = Mathieu(phi, morder, ellip, 'se');
  case 'even'
    A = Mathieu(phi, morder, ellip, 'ce');
  otherwise
    error('Unknown parity value, must be even or odd');
end

% Calculate weights
lmode = -Lmax:Lmax;
bval = exp(1i.*lmode.*phi.');
weights = bval \ A.';

%% Calculate Bessel basis functions

% Use the BscBessel function from OTTv1, requires OTT to be on the path.
basis = ott.BscBessel(Nmax, theta, 'lmode', lmode);

% In OTTv2:
% beam = basis .* weights.';

% In OTTv1:
beam = basis * weights;

%% Generate a visualisation of the beam

rng = [5, 5];   % Range (units:  medium wavelengths)

figure();
beam.visualise('range', rng);
colormap hot;
