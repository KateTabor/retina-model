%% Inspect spectra of side vs center bar in color-bar image
clear; close all;

% --- 1. Path to XYZ float TIFF ---
tiffPath = '/Users/kate/Documents/retina-model/bars_by1_b+100_XYZ.tiff';

% --- 2. Build scene from XYZ TIFF (uses helper function) ---
sceneFOVdeg = 2;   % same as before
scene = sceneFromMyXYZTiffv0(tiffPath, sceneFOVdeg);   % requires .m file on path

% --- 3. Get wavelength axis and photon cube from the scene ---
wave    = scene.spectrum.wave;      % [nWave x 1] wavelength samples
photons = scene.data.photons;       % [rows x cols x nWave]

[rows, cols, nWave] = size(photons);

% --- 4. Choose pixels: middle row, side bar, center bar ---
rowMid    = round(rows/2);   % vertical middle
colSide   = 50;              % side bar (50 px from left edge)
colCenter = 255;             % center bar

% Clamp column indices just in case
colSide   = max(1, min(cols, colSide));
colCenter = max(1, min(cols, colCenter));

% Extract spectra (photons vs wavelength) for those two pixels
specSide   = squeeze(photons(rowMid, colSide,   :));   % [nWave x 1]
specCenter = squeeze(photons(rowMid, colCenter, :));   % [nWave x 1]

% --- 5. Get an RGB image version of the scene for visualization ---
rgbImage = sceneGet(scene, 'rgb image');   % ISETBio helper

% --- 6. Plot: image on the left, spectra on the right ---
figure;

% Left: original color-bar image
subplot(1,2,1);
image(rgbImage);
axis image off;
title('Original color-bar image');

% Right: photon spectra for side vs center bar
subplot(1,2,2);
plot(wave, specSide,   '-o'); hold on;
plot(wave, specCenter, '-o');
xlabel('Wavelength (nm)');
ylabel('Photons (arb units)');
legend( ...
    sprintf('Side bar (row %d, col %d)',   rowMid, colSide), ...
    sprintf('Center bar (row %d, col %d)', rowMid, colCenter), ...
    'Location', 'best');
title('Spectra at side vs center bar');
grid on;
