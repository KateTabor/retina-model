%% Cone-surround *counts* visualization for LM vs LMS  (thresholded)

% --- Connectivity matrices
C_center_LM   = mLM.rgcRFcenterConeConnectivityMatrix;   
S_surround_LM = mLM.rgcRFsurroundConeConnectivityMatrix;

C_center_LMS   = mLMS.rgcRFcenterConeConnectivityMatrix;  
S_surround_LMS = mLMS.rgcRFsurroundConeConnectivityMatrix;

% --- Cone mosaics and indices
cm_LM  = mLM.inputConeMosaic;
cm_LMS = mLMS.inputConeMosaic;

lIdx_LM = cm_LM.lConeIndices;
mIdx_LM = cm_LM.mConeIndices;
sIdx_LM = cm_LM.sConeIndices;

lIdx_LMS = cm_LMS.lConeIndices;
mIdx_LMS = cm_LMS.mConeIndices;
sIdx_LMS = cm_LMS.sConeIndices;

% --- Surround weights restricted to each cone class
W_S_LM  = full(S_surround_LM(sIdx_LM,:));   % S-cone rows
W_M_LM  = full(S_surround_LM(mIdx_LM,:));   % M-cone rows
W_L_LM  = full(S_surround_LM(lIdx_LM,:));   % L-cone rows

W_S_LMS = full(S_surround_LMS(sIdx_LMS,:));
W_M_LMS = full(S_surround_LMS(mIdx_LMS,:));
W_L_LMS = full(S_surround_LMS(lIdx_LMS,:));

% --- Use each mosaic's *own* surround-threshold
thr_LM  = mLM.minSurroundWeightForInclusionInComputing;
thr_LMS = mLMS.minSurroundWeightForInclusionInComputing;

% --- COUNT of cones (above threshold) per RGC
numS_LM  = sum(abs(W_S_LM)  >= thr_LM, 1);
numM_LM  = sum(abs(W_M_LM)  >= thr_LM, 1);
numL_LM  = sum(abs(W_L_LM)  >= thr_LM, 1);

numS_LMS = sum(abs(W_S_LMS) >= thr_LMS, 1);
numM_LMS = sum(abs(W_M_LMS) >= thr_LMS, 1);
numL_LMS = sum(abs(W_L_LMS) >= thr_LMS, 1);

% --- RGC positions
pos_LM  = mLM.rgcRFpositionsDegs;   % [nRGC x 2]
pos_LMS = mLMS.rgcRFpositionsDegs;

% --- Shared color scale across ALL subplots
allCounts = [numS_LM numM_LM numL_LM numS_LMS numM_LMS numL_LMS];
cmin = 0;
cmax = max(allCounts);
if cmax == 0
    cmax = 1;   % avoid degenerate color axis
end

% --- Shared axis limits across panels (LM+LMS)
xAll = [pos_LM(:,1); pos_LMS(:,1)];
yAll = [pos_LM(:,2); pos_LMS(:,2)];
xlims = [min(xAll) max(xAll)];
ylims = [min(yAll) max(yAll)];

% --- Figure: 3 rows (S/M/L), 2 cols (LM / LMS)
figure;
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
colormap(parula);

% ---- Row 1: S-cone surround counts
nexttile;
scatter(pos_LM(:,1), pos_LM(:,2), 15, numS_LM, 'filled');
axis equal; xlim(xlims); ylim(ylims);
set(gca,'CLim',[cmin cmax]);
title('LM: S surround count');
xlabel('deg'); ylabel('deg');

nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, numS_LMS, 'filled');
axis equal; xlim(xlims); ylim(ylims);
set(gca,'CLim',[cmin cmax]);
title('LMS: S surround count');
xlabel('deg'); ylabel('deg');

% ---- Row 2: M-cone surround counts
nexttile;
scatter(pos_LM(:,1), pos_LM(:,2), 15, numM_LM, 'filled');
axis equal; xlim(xlims); ylim(ylims);
set(gca,'CLim',[cmin cmax]);
title('LM: M surround count');
xlabel('deg'); ylabel('deg');

nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, numM_LMS, 'filled');
axis equal; xlim(xlims); ylim(ylims);
set(gca,'CLim',[cmin cmax]);
title('LMS: M surround count');
xlabel('deg'); ylabel('deg');

% ---- Row 3: L-cone surround counts
nexttile;
scatter(pos_LM(:,1), pos_LM(:,2), 15, numL_LM, 'filled');
axis equal; xlim(xlims); ylim(ylims);
set(gca,'CLim',[cmin cmax]);
title('LM: L surround count');
xlabel('deg'); ylabel('deg');

nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, numL_LMS, 'filled');
axis equal; xlim(xlims); ylim(ylims);
set(gca,'CLim',[cmin cmax]);
title('LMS: L surround count');
xlabel('deg'); ylabel('deg');

% --- One shared colorbar
cb = colorbar;
cb.Layout.Tile = 'east';     % put it to the right of all tiles
cb.Label.String = 'Number of cones in surround (|w| ≥ threshold)';
