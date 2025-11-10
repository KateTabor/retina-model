function loop_enableS_by_coords
% loop_enableS_by_coords
% Run 4 conditions (0, 5, 50, 100% S-surround enabling) on the
% 2-deg prebaked ON mRGC mosaic, using the authors' own
% "give me coordinates → I will find nearest RGCs" style.
%
% - keeps the authors' console printout
% - keeps mRGC noise (we do NOT set noiseFlag = 'none')
% - after each condition, collapses the response to ONE number per RGC
%   (mean over time) **preferring the NOISY response**
% - saves everything to mRGC_apple_responses.mat (-v7.3)
% - makes a grouped bar plot for 5 cells

ieInit;

%% 1) build the stimulus / OI 
tiffPath = '/Users/kate/Documents/retina-model/image-set/train/trial_000001/initial_XYZ.tiff';
scene    = sceneFromMyXYZTiffv0(tiffPath, 2);
scene    = sceneInterpolate(scene, [256 256]);
oi       = oiCompute(oiCreate('wvf human'), scene);

%% 2) list of S-surround conditions we want
percentages = [0, 0.05, 0.50, 1.00];

%% 3) load ON mRGC mosaic ONCE to measure footprint
coneMosaicSpecies    = 'human';
opticsSubjectName    = 'PLOSpaperDefaultSubject';
rgcMosaicName        = 'PLOSpaperTemporal2DegsMosaic';
targetVisualSTFdescr = 'default';

baseM = mRGCMosaic.loadPrebakedMosaic( ...
    coneMosaicSpecies, ...
    opticsSubjectName, ...
    rgcMosaicName, ...
    targetVisualSTFdescr, ...
    'computeTheMosaicOptics', false);

allPos = baseM.rgcRFpositionsDegs;   % [N x 2]
nRGCs  = baseM.rgcsNum;

minX = min(allPos(:,1));  maxX = max(allPos(:,1));
minY = min(allPos(:,2));  maxY = max(allPos(:,2));
widthDeg  = maxX - minX;
heightDeg = maxY - minY;
centerX   = mean(allPos(:,1));
centerY   = mean(allPos(:,2));

fprintf('Base mosaic: N=%d, X=[%.2f %.2f], Y=[%.2f %.2f]\n', ...
    nRGCs, minX, maxX, minY, maxY);

% make the disk
radiusDeg = 0.9 * min(widthDeg, heightDeg);

% reproducible coordinates
rng(17);

% where we stash per-condition output
results = struct([]);

for ii = 1:numel(percentages)
    pct = percentages(ii);

    % ------------------------------------------------------------
    % A. load a CLEAN mosaic for THIS condition
    %     (this keeps the authors' "Loading prebaked..." print)
    % ------------------------------------------------------------
    theMRGCmosaic = mRGCMosaic.loadPrebakedMosaic( ...
        coneMosaicSpecies, ...
        opticsSubjectName, ...
        rgcMosaicName, ...
        targetVisualSTFdescr, ...
        'computeTheMosaicOptics', false);

    thisPos = theMRGCmosaic.rgcRFpositionsDegs;
    thisN   = theMRGCmosaic.rgcsNum;

    % ------------------------------------------------------------
    % B. choose which RGCs to enable
    % ------------------------------------------------------------
    if pct == 0
        % baseline: no changes
        rgcIdx     = [];
        nAttempted = 0;
        msg        = '';
    elseif pct >= 0.999
        % 100% condition: deterministic, all cells
        callStr = [ ...
            'theMRGCmosaic.enableSconeSurroundPoolingInSelectCells(' ...
            '1:' num2str(thisN) ...
            ', ''visualizeBeforeAndAfterConePoolingMaps'', false, ' ...
            '''exportVisualizationPDFdirectory'', '''');' ...
            ];
        msg = evalc(callStr);
        rgcIdx     = 1:thisN;
        nAttempted = thisN;
    else
        % 5% or 50%: authors' style — make random coords in a big disk
        nTarget = round(pct * thisN);

        theta  = 2*pi*rand(nTarget,1);
        r      = radiusDeg * sqrt(rand(nTarget,1));
        coords = [centerX + r.*cos(theta), ...
                  centerY + r.*sin(theta)];

        % map coords → nearest RGCs
        [~, rgcIdx] = pdist2(thisPos, coords, 'euclidean', 'Smallest', 1);
        rgcIdx = unique(rgcIdx);                      % remove dups
        nAttempted = numel(rgcIdx);

        callStr = sprintf([ ...
            'theMRGCmosaic.enableSconeSurroundPoolingInSelectCells(%s, ', ...
            '''visualizeBeforeAndAfterConePoolingMaps'', false, ', ...
            '''exportVisualizationPDFdirectory'', '''');' ...
            ], mat2str(rgcIdx));
        msg = evalc(callStr);
    end

    % ------------------------------------------------------------
    % C. classify what happened using the authors' own messages
    % ------------------------------------------------------------
    nEnabled = countSubstring(msg, 'Enabling S-cone inputs to RF surround of');
    nNoS     = countSubstring(msg, 'No S-cones in the input cone mosaic surrounding mRGC');
    nOther   = max(0, nAttempted - (nEnabled + nNoS));

    % ------------------------------------------------------------
    % D. compute THIS mosaic's response to THIS OI
    %     (keep NOISE; do not touch noiseFlag)
    % ------------------------------------------------------------
    % find its internal cone mosaic
    cm_for_rgc = [];
    for propName = {'inputConeMosaic','coneMosaic','theConeMosaic','theInputConeMosaic'}
        if isprop(theMRGCmosaic, propName{1})
            cm_for_rgc = theMRGCmosaic.(propName{1});
            break;
        end
    end
    if isempty(cm_for_rgc)
        error('No internal cone mosaic found for this mRGC mosaic.');
    end

    cmResp_internal = cm_for_rgc.compute(oi);
    cmResp_internal = double(squeeze(cmResp_internal));

    % normalize cone response into [1 x T x C]
    if isvector(cmResp_internal)
        nCones = numel(cmResp_internal);
        coneResp3D = reshape(cmResp_internal, [1 1 nCones]);
    elseif ndims(cmResp_internal) == 2
        [d1,d2] = size(cmResp_internal);
        if d1 > 1 && d2 > 1
            nCones = d1; nTime = d2;
            coneResp3D = reshape(cmResp_internal.', [1 nTime nCones]);
        else
            nCones = max(d1,d2);
            coneResp3D = reshape(cmResp_internal, [1 1 nCones]);
        end
    else
        [rS,cS,nTime] = size(cmResp_internal);
        nCones = rS * cS;
        tmp = reshape(cmResp_internal, [nCones, nTime]);
        coneResp3D = reshape(tmp.', [1 nTime nCones]);
    end

    % tile in time if single frame
    if size(coneResp3D,2) == 1
        coneResp3D = repmat(coneResp3D, [1 50 1]);
    end

    % time support
    if isprop(cm_for_rgc, 'integrationTime') && ~isempty(cm_for_rgc.integrationTime)
        dt = cm_for_rgc.integrationTime;
    else
        dt = 0.010;
    end
    coneRespTime = (0:(size(coneResp3D,2)-1)) * dt;

    % actual mRGC compute (WITH noise)
    [rgcResp, rgcRespNoisy, rgcTime] = theMRGCmosaic.compute( ...
        coneResp3D, ...
        coneRespTime);

    % ------------------------------------------------------------
    % E. collapse over time → ONE number per cell
    %     *prefer the NOISY response*
    % ------------------------------------------------------------
    if ~isempty(rgcRespNoisy)
        % rgcRespNoisy: [nTrials x T x N]
        respToUse = squeeze(rgcRespNoisy);  % → [T x N] if nTrials=1
    else
        % fall back to noise-free
        respToUse = squeeze(rgcResp);       % → [T x N]
    end

    if isvector(respToUse)
        % weird edge case
        perCellMean = respToUse(:).';       % 1 x N
    else
        perCellMean = mean(respToUse, 1);   % 1 x N
    end

    % ------------------------------------------------------------
    % F. stash everything
    % ------------------------------------------------------------
    results(ii).pct         = pct;
    results(ii).nRGCs       = thisN;
    results(ii).nAttempted  = nAttempted;
    results(ii).nEnabled    = nEnabled;
    results(ii).nNoS        = nNoS;
    results(ii).nOther      = nOther;
    results(ii).rgcIdxUsed  = rgcIdx;
    results(ii).rgcResp     = rgcResp;
    results(ii).rgcRespNoisy= rgcRespNoisy;
    results(ii).rgcTime     = rgcTime;
    results(ii).perCellMean = perCellMean;
end

%% 4) short console summary
for ii = 1:numel(results)
    r = results(ii);
    fprintf('pct=%3.0f%%  attempted=%5d  enabled=%5d  no-S=%5d  other=%5d  (N=%d)\n', ...
        100*r.pct, r.nAttempted, r.nEnabled, r.nNoS, r.nOther, r.nRGCs);
end

%% 5) save to MAT (for the AI model)
save('mRGC_apple_responses.mat', 'percentages', 'results', '-v7.3');
fprintf('Saved all responses to mRGC_apple_responses.mat\n');

%% 6) BAR PLOT for 5 cells
% pick cells that 5% actually enabled
enabled5 = results(2).rgcIdxUsed;

if isempty(enabled5)
    cellsToPlot = 1:5;  % fallback
    fprintf('5%% enabled 0 cells; plotting cells 1:5 instead.\n');
else
    nTake = min(5, numel(enabled5));
    cellsToPlot = enabled5(1:nTake);
    fprintf('Plotting cells (5%%-enabled): %s\n', mat2str(cellsToPlot));
end

nConds = numel(results);
nCells = numel(cellsToPlot);
vals   = nan(nCells, nConds);   % rows = cells, cols = conditions

for ci = 1:nCells
    cID = cellsToPlot(ci);
    for jj = 1:nConds
        r = results(jj);
        if cID <= numel(r.perCellMean)
            vals(ci, jj) = r.perCellMean(cID);
        else
            vals(ci, jj) = NaN;
        end
    end
end

figure; clf;
bar(1:nCells, vals, 'grouped');
xticks(1:nCells);
xticklabels(arrayfun(@(x) sprintf('RGC %d', x), cellsToPlot, 'UniformOutput', false));
xlabel('Cells');
ylabel('mean mRGC response (noisy, averaged over time)');
legend({'0%','5%','50%','100%'}, 'Location', 'bestoutside');
title('mRGC responses across S-surround conditions (5 cells)');
grid on;

end  % main function


% ------------------------------------------------------------
% helper: count occurrences of a substring
% ------------------------------------------------------------
function n = countSubstring(str, pattern)
if isempty(str)
    n = 0;
else
    n = numel(strfind(str, pattern));
end
end
