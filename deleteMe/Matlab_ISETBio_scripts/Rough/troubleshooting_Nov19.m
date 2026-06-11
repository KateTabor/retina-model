% Make sure cm exists (e.g., from build2degONmRGCMosaic_v2)
figure; hold on;

% Plot L, M, S cone quantal efficiencies vs wavelength
plot(cm.wave, cm.qe(:, cm.LCONE_ID), 'r-',  'LineWidth', 2); % L (red)
plot(cm.wave, cm.qe(:, cm.MCONE_ID), 'g-',  'LineWidth', 2); % M (green)
plot(cm.wave, cm.qe(:, cm.SCONE_ID), 'b-',  'LineWidth', 2); % S (blue)

xlabel('Wavelength (nm)');
ylabel('Quantal efficiency (arb. units)');
title('ISETBio cone spectral sensitivities');
legend({'L cones','M cones','S cones'}, 'Location', 'best');
grid on;
hold off;

%%%%%%%%%
%% Simple cone response visualization (S, M, L) with equal-sized plots

% Get cone mosaic from the LM mRGC mosaic
cm = mLM.inputConeMosaic;

% Cone positions (deg) and indices by type
pos  = cm.coneRFpositionsDegs;   % [Ncones x 2]
Lidx = cm.lConeIndices;
Midx = cm.mConeIndices;
Sidx = cm.sConeIndices;

% Cone responses from runOneMosaic (single frame)
resp = coneRespSingle(:);        % make sure it's a column

% Shared color limits across S/M/L plots
vmin = min(resp);
vmax = max(resp);

% Get an RGB image representation of the scene
rgb = sceneGet(scene, 'rgb image');

figure;
t = tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact');

%% 1. Original image
ax1 = nexttile;
image(rgb);
axis image off;
title('Original image');

%% 2. S-cone responses
ax2 = nexttile;
scatter(pos(Sidx,1), pos(Sidx,2), 10, resp(Sidx), 'filled');
axis equal tight;
set(gca,'YDir','reverse');
xlabel('deg'); ylabel('deg');
title('S-cone responses');
clim([vmin vmax]);

%% 3. M-cone responses
ax3 = nexttile;
scatter(pos(Midx,1), pos(Midx,2), 10, resp(Midx), 'filled');
axis equal tight;
set(gca,'YDir','reverse');
xlabel('deg'); ylabel('deg');
title('M-cone responses');
clim([vmin vmax]);

%% 4. L-cone responses
ax4 = nexttile;
scatter(pos(Lidx,1), pos(Lidx,2), 10, resp(Lidx), 'filled');
axis equal tight;
set(gca,'YDir','reverse');
xlabel('deg'); ylabel('deg');
title('L-cone responses');
clim([vmin vmax]);

% Shared colorbar for the tiled layout (keeps axes the same size)
cb = colorbar(ax4, 'eastoutside');
ylabel(cb, 'Cone response (a.u.)');

colormap(parula);   % one colormap for S/M/L plots
