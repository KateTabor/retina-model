%% plot_bars_allFigures.m
% Visualize all bars_gray1_* images and their LM / LMS / LMS_norm responses
% -- or bars_rg1_*

clear; close all;

%% 1. Directory with response .mat files
respDir = '/Users/kate/Documents/retina-model/2025-11-18_bars_gradients/20251124_runAllImages';

% Find the five bars_rg1_* response files
respFiles = dir(fullfile(respDir, 'mRGCresp_*bars_rg1_*_XYZ_dt10ms_T10.mat'));
if numel(respFiles) ~= 5
    warning('Expected 5 bars_rg1_ files; found %d', numel(respFiles));
end

% Sort by name so we get a+000, a+050, ..., a-100 in order
[~, idx]  = sort({respFiles.name});
respFiles = respFiles(idx);

%% 2. Load responses, image paths, positions, scenes, and global ranges

data    = struct([]);
allVals = [];
allPos  = [];

for ii = 1:numel(respFiles)
    D = load(fullfile(respDir, respFiles(ii).name));

    % Store responses as column vectors
    data(ii).respLM        = D.respLM(:);
    data(ii).respLMS       = D.respLMS(:);
    data(ii).respLMS_norm  = D.respLMS_norm(:);
    data(ii).thisImg       = D.thisImg;
    data(ii).rgc_positions = D.rgc_positions;

    % Build and store scene + image name (for titles)
    sceneFOVdeg     = 2;   % same as used when generating responses
    data(ii).scene  = sceneFromXYZfloat32(data(ii).thisImg, sceneFOVdeg);
    [~, imgName, ext] = fileparts(data(ii).thisImg);
    data(ii).imgName = [imgName ext];

    % Accumulate values for global color limits / line-plot limits
    allVals = [allVals; ...
               data(ii).respLM; ...
               data(ii).respLMS; ...
               data(ii).respLMS_norm]; 

    % Positions (same mosaic across files, but safe to collect)
    allPos = [allPos; data(ii).rgc_positions]; 
end

% Number of rows/images to plot (for Figures 1–3)
nRows = numel(data);

% Global histogram bins and color limits (used for maps / line ylim)
nBins = 60;
edges = linspace(min(allVals), max(allVals), nBins); %#ok<NASGU>
cmin  = min(allVals); %#ok<NASGU>
cmax  = max(allVals); %#ok<NASGU>

respMin = min(allVals);
respMax = max(allVals);

% Global x/y limits for scatter maps (Figure 2)
xlims = [min(allPos(:,1)) max(allPos(:,1))];
ylims = [min(allPos(:,2)) max(allPos(:,2))];

% Colors
LMcolor      = [0.8 0.1 0.1];  % red
LMScolor     = [0.1 0.1 0.8];  % blue
LMSnormColor = [0.0 0.6 0.6];  % cyan

% Tolerance for picking "mid-row" in Figure 3 (in deg)
tolLine = 0.02;

%% FIGURE 1: nRows×3 (image, histogram, scatter)

figure(1); clf;
t1 = tiledlayout(nRows, 3, ...
    'TileSpacing','tight', ...
    'Padding','tight');

for rr = 1:nRows
    d = data(rr);

    % Row-specific min/max across LM, LMS, LMS_norm responses
    rowVals = [d.respLM; d.respLMS; d.respLMS_norm];
    rowMin  = min(rowVals);
    rowMax  = max(rowVals);
    rowEdges = linspace(rowMin, rowMax, 50);

    % --- Column 1: original image
    nexttile(t1);
    rgb = sceneGet(d.scene, 'rgb image');
    image(rgb);
    axis image off;
    title(d.imgName, 'Interpreter','none');

    % --- Column 2: overlaid histograms (LM, LMS, LMS_norm)
    axHist = nexttile(t1);
    hold(axHist, 'on');
    histogram(axHist, d.respLM,       rowEdges, ...
              'FaceColor', LMcolor, ...
              'FaceAlpha', 0.4, ...
              'EdgeColor', 'none');
    histogram(axHist, d.respLMS,      rowEdges, ...
              'FaceColor', LMScolor, ...
              'FaceAlpha', 0.4, ...
              'EdgeColor', 'none');
    histogram(axHist, d.respLMS_norm, rowEdges, ...
              'FaceColor', LMSnormColor, ...
              'FaceAlpha', 0.4, ...
              'EdgeColor', 'none');
    hold(axHist, 'off');
    xlabel('Response (a.u.)');
    ylabel('Count');
    title('RGC response hist');
    xlim(axHist, [rowMin rowMax]);   % per-row x-range

    % Legend only on the top row histogram
    if rr == 1
        legend(axHist, {'LM','LMS','LMS norm'}, 'Location','best');
    end

    % --- Column 3: scatter LMS vs LMS_norm with matched axes
    axScat = nexttile(t1);
    scatter(axScat, d.respLMS, d.respLMS_norm, 5, 'filled');
    xlabel('respLMS');
    ylabel('respLMS\_norm');
    title('LMS vs LMS\_norm per RGC');
    grid(axScat, 'on');

    % Make scatter x- and y-axis match the histogram x-axis range
    xlim(axScat, [rowMin rowMax]);
    ylim(axScat, [rowMin rowMax]);
    axis(axScat, 'square');
end

% Make all axes text slightly smaller in this figure
ax1 = findall(gcf, 'Type','axes');
set(ax1, 'FontSize', 8);

%% FIGURE 2: nRows×4 maps (image, LM, LMS, LMS_norm)

figure(2); clf;
t2 = tiledlayout(nRows, 4, ...
    'TileSpacing','tight', ...
    'Padding','tight');

for rr = 1:nRows
    d = data(rr);

    % Column 1: image
    nexttile(t2);
    rgb = sceneGet(d.scene, 'rgb image');
    image(rgb);
    axis image off;
    title(d.imgName, 'Interpreter','none');

    % Column 2: LM map
    axLM = nexttile(t2);
    scatter(axLM, d.rgc_positions(:,1), d.rgc_positions(:,2), 5, d.respLM, 'filled');
    axis(axLM, 'equal');
    xlim(axLM, xlims); ylim(axLM, ylims);
    title('LM');
    clim(axLM, [respMin respMax]);
    xlabel('deg'); ylabel('deg');

    % Column 3: LMS map
    axLMS = nexttile(t2);
    scatter(axLMS, d.rgc_positions(:,1), d.rgc_positions(:,2), 5, d.respLMS, 'filled');
    axis(axLMS, 'equal');
    xlim(axLMS, xlims); ylim(axLMS, ylims);
    title('LMS');
    clim(axLMS, [respMin respMax]);
    xlabel('deg'); ylabel('deg');

    % Column 4: LMS_norm map
    axNorm = nexttile(t2);
    scatter(axNorm, d.rgc_positions(:,1), d.rgc_positions(:,2), 5, d.respLMS_norm, 'filled');
    axis(axNorm, 'equal');
    xlim(axNorm, xlims); ylim(axNorm, ylims);
    title('LMS\_norm');
    clim(axNorm, [respMin respMax]);
    xlabel('deg'); ylabel('deg');
end

colormap(parula);

% Colorbar for all response maps
cb2 = colorbar;
cb2.Layout.Tile = 'east';
ylabel(cb2,'Response (a.u.)');

% Make all axes text slightly smaller in this figure
ax2 = findall(gcf, 'Type','axes');
set(ax2, 'FontSize', 8);

%% FIGURE 3: nRows×2 (image, line plot at mid-row)

figure(3); clf;
t3 = tiledlayout(nRows, 2, ...
    'TileSpacing','tight', ...
    'Padding','tight');

for rr = 1:nRows
    d = data(rr);

    % --- Column 1: image
    nexttile(t3);
    rgb = sceneGet(d.scene, 'rgb image');
    image(rgb);
    axis image off;
    title(d.imgName, 'Interpreter','none');

    % --- Column 2: line plot of LM/LMS/LMS_norm vs x-position at mid-row
    axLine = nexttile(t3); 
    hold(axLine, 'on');

    % Compute mid-row indices for THIS dataset
    yPos = d.rgc_positions(:,2);
    yMid = median(yPos);

    midIdx = find(abs(yPos - yMid) < tolLine);
    if isempty(midIdx)
        % Fallback: take the single closest row
        [~, closestIdx] = min(abs(yPos - yMid));
        midIdx = closestIdx;
    end

    % Sort by x-position along that mid-row
    xRow = d.rgc_positions(midIdx, 1);
    [xSorted, order] = sort(xRow);

    lm_line       = d.respLM(midIdx(order));
    lms_line      = d.respLMS(midIdx(order));
    lms_norm_line = d.respLMS_norm(midIdx(order));

    % Plot
    hLM       = plot(axLine, xSorted, lm_line,       'Color', LMcolor,      'LineWidth', 1.0);
    hLMS      = plot(axLine, xSorted, lms_line,      'Color', LMScolor,     'LineWidth', 1.0);
    hLMS_norm = plot(axLine, xSorted, lms_norm_line, 'Color', LMSnormColor, 'LineWidth', 1.0);

    xlabel('x-position (deg)');
    ylabel('RGC response (a.u.)');
    title('Mid-row RGC responses');
    grid(axLine, 'on');

    % Legend only on the top row line plot
    if rr == 1
        legend(axLine, [hLM, hLMS, hLMS_norm], ...
            {'LM','LMS','LMS norm'}, 'Location','best');
    end

    % Per-row vertical range based on this row's mid-line values
    midVals = [lm_line(:); lms_line(:); lms_norm_line(:)];
    yMin = min(midVals);
    yMax = max(midVals);
    pad = 0.05 * (yMax - yMin + eps);
    ylim(axLine, [yMin - pad, yMax + pad]);

    hold(axLine, 'off');
end

% Make all axes text slightly smaller in this figure
ax3 = findall(gcf, 'Type','axes');
set(ax3, 'FontSize', 8);



% 
% %% plot_bars_allFigures.m
% % Visualize all bars_gray1_* images and their LM / LMS / LMS_norm responses
% % -- or bars_by1_
% 
% clear; % close all;
% 
% %% 1. Directory with response .mat files
% respDir = '/Users/kate/Documents/retina-model/2025-11-18_bars_gradients/20251124_runAllImages';
% 
% % Find the five bars_rg1_* response files
% respFiles = dir(fullfile(respDir, 'mRGCresp_*bars_rg1_*_XYZ_dt10ms_T10.mat'));
% if numel(respFiles) ~= 5
%     warning('Expected 5 bars_rg1_ files; found %d', numel(respFiles));
% end
% 
% % Sort by name so we get L000, L025, ..., L100 in order
% [~, idx]  = sort({respFiles.name});
% respFiles = respFiles(idx);
% 
% %% 2. Load responses, image paths, positions, scenes, and global ranges
% 
% data    = struct([]);
% allVals = [];
% allPos  = [];
% 
% for ii = 1:numel(respFiles)
%     D = load(fullfile(respDir, respFiles(ii).name));
% 
%     % Store responses as column vectors
%     data(ii).respLM        = D.respLM(:);
%     data(ii).respLMS       = D.respLMS(:);
%     data(ii).respLMS_norm  = D.respLMS_norm(:);
%     data(ii).thisImg       = D.thisImg;
%     data(ii).rgc_positions = D.rgc_positions;
% 
%     % Build and store scene + image name (for titles)
%     sceneFOVdeg     = 2;   % same as used when generating responses
%     data(ii).scene  = sceneFromXYZfloat32(data(ii).thisImg, sceneFOVdeg);
%     [~, imgName, ext] = fileparts(data(ii).thisImg);
%     data(ii).imgName = [imgName ext];
% 
%     % Accumulate values for global color limits / line-plot limits
%     allVals = [allVals; ...
%                data(ii).respLM; ...
%                data(ii).respLMS; ...
%                data(ii).respLMS_norm]; 
% 
%     % Positions (same mosaic across files, but safe to collect)
%     allPos = [allPos; data(ii).rgc_positions]; 
% 
%     % Precompute mid-row line data for Figure 3
%     yPos = data(ii).rgc_positions(:,2);
%     yMid = median(yPos);
%     tol  = 0.02;   % deg tolerance around mid-row
% 
%     midIdx = find(abs(yPos - yMid) < tol);
%     if isempty(midIdx)
%         % Fallback: take the single closest row
%         [~, closestIdx] = min(abs(yPos - yMid));
%         midIdx = closestIdx;
%     end
% 
%     [xSorted, order] = sort(data(ii).rgc_positions(midIdx,1));
%     data(ii).xPos              = xSorted;
%     data(ii).respLM_line       = data(ii).respLM(midIdx(order));
%     data(ii).respLMS_line      = data(ii).respLMS(midIdx(order));
%     data(ii).respLMS_norm_line = data(ii).respLMS_norm(midIdx(order));
% end
% 
% % Number of rows/images to plot (for Figures 1–3)
% nRows = numel(data);
% 
% % Global histogram bins and color limits (used for maps / line ylim)
% nBins = 60;
% edges = linspace(min(allVals), max(allVals), nBins); 
% cmin  = min(allVals); 
% cmax  = max(allVals); 
% 
% respMin = min(allVals);
% respMax = max(allVals);
% 
% % Global x/y limits for scatter maps (Figure 2)
% xlims = [min(allPos(:,1)) max(allPos(:,1))];
% ylims = [min(allPos(:,2)) max(allPos(:,2))];
% 
% % Colors
% LMcolor      = [0.8 0.1 0.1];  % red
% LMScolor     = [0.1 0.1 0.8];  % blue
% LMSnormColor = [0.0 0.6 0.6];  % cyan
% 
% %% FIGURE 1: nRows×3 (image, histogram, scatter)
% 
% figure(1); clf;
% t1 = tiledlayout(nRows, 3, ...
%     'TileSpacing','tight', ...
%     'Padding','tight');
% 
% for rr = 1:nRows
%     d = data(rr);
% 
%     % Row-specific min/max across LM, LMS, LMS_norm responses
%     rowVals = [d.respLM; d.respLMS; d.respLMS_norm];
%     rowMin  = min(rowVals);
%     rowMax  = max(rowVals);
%     rowEdges = linspace(rowMin, rowMax, 50);
% 
%     % --- Column 1: original image
%     nexttile(t1);
%     rgb = sceneGet(d.scene, 'rgb image');
%     image(rgb);
%     axis image off;
%     title(d.imgName, 'Interpreter','none');
% 
%     % --- Column 2: overlaid histograms (LM, LMS, LMS_norm)
%     axHist = nexttile(t1);
%     hold(axHist, 'on');
%     histogram(axHist, d.respLM,       rowEdges, ...
%               'FaceColor', LMcolor, ...
%               'FaceAlpha', 0.4, ...
%               'EdgeColor', 'none');
%     histogram(axHist, d.respLMS,      rowEdges, ...
%               'FaceColor', LMScolor, ...
%               'FaceAlpha', 0.4, ...
%               'EdgeColor', 'none');
%     histogram(axHist, d.respLMS_norm, rowEdges, ...
%               'FaceColor', LMSnormColor, ...
%               'FaceAlpha', 0.4, ...
%               'EdgeColor', 'none');
%     hold(axHist, 'off');
%     xlabel('Response (a.u.)');
%     ylabel('Count');
%     title('RGC response hist');
%     xlim(axHist, [rowMin rowMax]);   % per-row x-range
% 
%     % Legend only on the top row histogram
%     if rr == 1
%         legend(axHist, {'LM','LMS','LMS norm'}, 'Location','best');
%     end
% 
%     % --- Column 3: scatter LMS vs LMS_norm with matched axes
%     axScat = nexttile(t1);
%     scatter(axScat, d.respLMS, d.respLMS_norm, 5, 'filled');
%     xlabel('respLMS');
%     ylabel('respLMS\_norm');
%     title('LMS vs LMS\_norm per RGC');
%     grid(axScat, 'on');
% 
%     % Make scatter x- and y-axis match the histogram x-axis range
%     xlim(axScat, [rowMin rowMax]);
%     ylim(axScat, [rowMin rowMax]);
%     axis(axScat, 'square');
% end
% 
% % Make all axes text slightly smaller in this figure
% ax1 = findall(gcf, 'Type','axes');
% set(ax1, 'FontSize', 8);
% 
% %% FIGURE 2: nRows×4 maps (image, LM, LMS, LMS_norm)
% 
% figure(2); clf;
% t2 = tiledlayout(nRows, 4, ...
%     'TileSpacing','tight', ...
%     'Padding','tight');
% 
% for rr = 1:nRows
%     d = data(rr);
% 
%     % Column 1: image
%     nexttile(t2);
%     rgb = sceneGet(d.scene, 'rgb image');
%     image(rgb);
%     axis image off;
%     title(d.imgName, 'Interpreter','none');
% 
%     % Column 2: LM map
%     axLM = nexttile(t2);
%     scatter(axLM, d.rgc_positions(:,1), d.rgc_positions(:,2), 5, d.respLM, 'filled');
%     axis(axLM, 'equal');
%     xlim(axLM, xlims); ylim(axLM, ylims);
%     title('LM');
%     clim(axLM, [respMin respMax]);
%     xlabel('deg'); ylabel('deg');
% 
%     % Column 3: LMS map
%     axLMS = nexttile(t2);
%     scatter(axLMS, d.rgc_positions(:,1), d.rgc_positions(:,2), 5, d.respLMS, 'filled');
%     axis(axLMS, 'equal');
%     xlim(axLMS, xlims); ylim(axLMS, ylims);
%     title('LMS');
%     clim(axLMS, [respMin respMax]);
%     xlabel('deg'); ylabel('deg');
% 
%     % Column 4: LMS_norm map
%     axNorm = nexttile(t2);
%     scatter(axNorm, d.rgc_positions(:,1), d.rgc_positions(:,2), 5, d.respLMS_norm, 'filled');
%     axis(axNorm, 'equal');
%     xlim(axNorm, xlims); ylim(axNorm, ylims);
%     title('LMS\_norm');
%     clim(axNorm, [respMin respMax]);
%     xlabel('deg'); ylabel('deg');
% end
% 
% colormap(parula);
% 
% % Colorbar for all response maps
% cb2 = colorbar;
% cb2.Layout.Tile = 'east';
% ylabel(cb2,'Response (a.u.)');
% 
% % Make all axes text slightly smaller in this figure
% ax2 = findall(gcf, 'Type','axes');
% set(ax2, 'FontSize', 8);
% 
% %% FIGURE 3: nRows×2 (image, line plot at mid-row)
% 
% figure(3); clf;
% t3 = tiledlayout(nRows, 2, ...
%     'TileSpacing','tight', ...
%     'Padding','tight');
% 
% for rr = 1:nRows
%     d = data(rr);
% 
%     % --- Column 1: image
%     nexttile(t3);
%     rgb = sceneGet(d.scene, 'rgb image');
%     image(rgb);
%     axis image off;
%     title(d.imgName, 'Interpreter','none');
% 
%     % --- Column 2: line plot of LM/LMS/LMS_norm vs x-position at mid-row
%     axLine = nexttile(t3); hold(axLine, 'on');
% 
%     hLM       = plot(axLine, d.xPos, d.respLM_line,       'Color', LMcolor,      'LineWidth', 1.0);
%     hLMS      = plot(axLine, d.xPos, d.respLMS_line,      'Color', LMScolor,     'LineWidth', 1.0);
%     hLMS_norm = plot(axLine, d.xPos, d.respLMS_norm_line, 'Color', LMSnormColor, 'LineWidth', 1.0);
% 
%     xlabel('x-position (deg)');
%     ylabel('RGC response (a.u.)');
%     title('Mid-row RGC responses');
%     grid(axLine, 'on');
% 
%     % Legend only on the top row line plot
%     if rr == 1
%         legend(axLine, [hLM, hLMS, hLMS_norm], ...
%             {'LM','LMS','LMS norm'}, 'Location','best');
%     end
% 
%     % Optional: keep same vertical range across all rows
%     ylim(axLine, [respMin respMax]);
% 
%     hold(axLine, 'off');
% end
% 
% % Make all axes text slightly smaller in this figure
% ax3 = findall(gcf, 'Type','axes');
% set(ax3, 'FontSize', 8);
