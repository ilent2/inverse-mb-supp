% example_paraxial_mb
%
% Generate paraxial Mathieu beam visualisations, similar to figure 2
% from
%
%   Inverse-Fourier Non-diffracting Beams for Optical Trapping
%   Martinez-Ruiz et al. [Journal TBD], 2020.
%
% Copyright Isaac Lenton, 2020

% Setup the figure
figure();

% Functions to help with visualisation
zm = 150;
zoom = @(im) im(zm:end-zm, zm:end-zm);
zm = 200;
zoom2 = @(im) im(zm:end-zm, zm:end-zm);

%% Generate paraxial Mathieu beam

ifocal = 11;  % Interfocal separation
scale = 50;   % Pattern scale (pixels)
q = 80;       % MB ellipticity
m = 5;        % MB m-mode

sigma = 1;    % Gaussian convolution (pixels)

Escalar = BscPmMathieu.calculate_scalar_field('even', ...
  m, scale, q, ifocal, 'Nres', 512);

subplot(1, 3, 1);
imagesc(zoom(abs(Escalar)));
axis image;
title('Mathieu Beam');

Efar = fftshift(fft2(Escalar));

subplot(1, 3, 2);
imagesc(zoom2(abs(Efar)));
axis image;
title('Frequency Space');

subplot(1, 3, 3);
imagesc(imgaussfilt(zoom2(abs(Efar)), sigma));
axis image;
title('Gaussian-filtered freq. space');

colormap hot;
