% example_inverse_mb.m - Demonstrate the BscPmMathieu class
%
% For more complete documentation on BscPmMathieu, checkout the
% documentation by typing `help BscPmMathieu` into the command window.
%
% Requires OTTv1.  Demonstrates the functionality used to generate
% most of the figures in
%
%   Inverse-Fourier Non-diffracting Beams for Optical Trapping
%   Martinez-Ruiz et al. [Journal TBD], 2020.
%
% Copyright Isaac Lenton, 2020

index_medium = 1;  % Refractive index in medium
NA = 1;           % Numerical aperture of lens
scale = 0.14;   % Scale relative to objective aperture

order = 3;         % Mathieu beam order (m)
parity = 'even';               % Mathieu beam parity (q)

ellipticity = 30;  % Mathieu beam helicity
ifocal = 0.1;      % Interfocal separation

Nmax = 30;
Nres = 256;

beam = BscPmMathieu(NA, order, parity, ...
    'ellipticity', ellipticity, 'scale', scale, ...
    'index_medium', index_medium, 'ifocal', ifocal, ...
    'Nres', Nres, 'Nmax', Nmax, 'projection', 'tanasin');
  
%% Generate visualisation

rng = [5, 5];   % Range in wavelengths.  Be careful setting this to large
                % values if Nmax is too small.  Might need Nmax=50 or >100.

figure();
beam.visualise('range', rng);
colormap hot;
  
  