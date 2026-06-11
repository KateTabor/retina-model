function [noiseFreeMRGCresponses, noisyMRGCresponseInstances, responseTemporalSupportSeconds] = ...
    compute_mRGC_withCustomSurround(mosaic, ...
                                    theInputConeMosaicResponse, ...
                                    theInputConeMosaicResponseTemporalSupportSeconds, ...
                                    surroundOverride, ...
                                    varargin)
% COMPUTE_MRGC_WITHCUSTOMSURROUND
%   Clone of the *linear* center–surround part of mRGCMosaic.compute, with
%   one extra feature: you can optionally override the surround wiring.
%
%   [noiseFree, noisy, tSupport] = compute_mRGC_withCustomSurround( ...
%       mosaic, coneResp, tSupport, surroundOverride, 'seed', 1)
%
%   Inputs
%       mosaic   : mRGCMosaic object (LM or LMS)
%       coneResp : [nTrials x nTime x nCones] cone responses
%       tSupport : [1 x nTime] or [nTime x 1] temporal support (seconds)
%       surroundOverride : [] to use mosaic.rgcRFsurroundConeConnectivityMatrix,
%                          or [nCones x nRGC] to override the surround wiring.
%
%   Outputs
%       noiseFreeMRGCresponses     : [nTrials x nTime x nRGC] noise-free RGC resp
%       noisyMRGCresponseInstances : [nTrials x nTime x nRGC] noisy RGC resp
%       responseTemporalSupportSeconds : copy of tSupport

    % ---------- 0. Parse optional args (we only care about 'seed') ----------
    p = inputParser;
    p.addParameter('seed', [], @isnumeric);
    p.parse(varargin{:});
    noiseSeed = p.Results.seed;

    % ---------- 1. Basic checks on cone response tensor ----------
    assert(ndims(theInputConeMosaicResponse) == 3, ...
        'coneResp must be [nTrials x nTime x nCones].');

    [nTrials, nTime, nCones] = size(theInputConeMosaicResponse);

    assert(nTime == numel(theInputConeMosaicResponseTemporalSupportSeconds), ...
        '2nd dim of coneResp (%d) must match length of tSupport (%d).', ...
        nTime, numel(theInputConeMosaicResponseTemporalSupportSeconds));

    responseTemporalSupportSeconds = theInputConeMosaicResponseTemporalSupportSeconds(:).'; %#ok<NASGU>

    % ---------- 2. Basic checks on mosaic wiring ----------
    S_mosaic = mosaic.rgcRFsurroundConeConnectivityMatrix;
    C_mosaic = mosaic.rgcRFcenterConeConnectivityMatrix;

    [nCones_center, nRGC_center]   = size(C_mosaic);
    [nCones_surround, nRGC_surround] = size(S_mosaic);

    assert(nCones_center == nCones_surround, ...
        'Center and surround matrices must have same # of cone rows.');
    assert(nRGC_center == nRGC_surround, ...
        'Center and surround matrices must have same # of RGC columns.');
    assert(nCones_center == nCones, ...
        ['# of cones in coneResp (size(:,:,3)=%d) must match # of rows in ' ...
         'center/surround matrices (%d).'], nCones, nCones_center);

    nRGC = mosaic.rgcsNum;
    assert(nRGC == nRGC_center, ...
        ['mosaic.rgcsNum (%d) must match # of columns in center/surround ' ...
         'matrices (%d).'], nRGC, nRGC_center);

    % ---------- 3. Surround override checks ----------
    useCustomSurround = ~isempty(surroundOverride);
    if useCustomSurround
        [nCones_override, nRGC_override] = size(surroundOverride);
        assert(nCones_override == nCones_center, ...
            ['surroundOverride must have %d rows (one per cone); got %d.'], ...
            nCones_center, nCones_override);
        assert(nRGC_override == nRGC_center, ...
            ['surroundOverride must have %d columns (one per RGC); got %d.'], ...
            nRGC_center, nRGC_override);
    end

    % ---------- 4. Preallocate output ----------
    noiseFreeMRGCresponses = zeros(nTrials, nTime, nRGC);

    % ---------- 5. Main center/surround computation ----------
    % (Use plain FOR instead of PARFOR to get clearer error messages while
    %  debugging. You can switch back to parfor once everything is solid.)
    for iRGC = 1:nRGC
        % ----- Center wiring -----
        centerConnectivityVector = full(C_mosaic(:, iRGC));
        centerConeIndices = find(centerConnectivityVector > ...
                                 mRGCMosaic.minCenterWeightForInclusionInComputing);
        centerConeWeights = reshape(centerConnectivityVector(centerConeIndices), ...
                                    [1 1 numel(centerConeIndices)]);

        % Spatial pooling over center cones
        centerSpatiallyIntegratedActivations = sum( ...
            bsxfun(@times, ...
                   theInputConeMosaicResponse(:,:,centerConeIndices), ...
                   centerConeWeights), ...
            3);

        % ----- Surround wiring -----
        if ~isempty(S_mosaic)
            if useCustomSurround
                surroundConnectivityVector = full(surroundOverride(:, iRGC));
            else
                surroundConnectivityVector = full(S_mosaic(:, iRGC));
            end

            surroundConeIndices = find(surroundConnectivityVector > ...
                mRGCMosaic.minSurroundWeightForInclusionInComputing);

            surroundConeWeights = reshape( ...
                surroundConnectivityVector(surroundConeIndices), ...
                [1 1 numel(surroundConeIndices)]);

            surroundSpatiallyIntegratedActivations = sum( ...
                bsxfun(@times, ...
                       theInputConeMosaicResponse(:,:,surroundConeIndices), ...
                       surroundConeWeights), ...
                3);
        else
            surroundSpatiallyIntegratedActivations = ...
                0 * centerSpatiallyIntegratedActivations;
        end

        % ----- Linear center - surround, scaled by responseGains -----
        noiseFreeMRGCresponses(:,:,iRGC) = mosaic.responseGains(iRGC) * ...
            (centerSpatiallyIntegratedActivations - ...
             surroundSpatiallyIntegratedActivations);
    end

    % ---------- 6. Noise handling (same logic as original compute) ----------
    if isempty(mosaic.noiseFlag)
        fprintf(['Warning: The mRGCMosaic.noiseFlag not set before calling ' ...
                 'compute_mRGC_withCustomSurround(). Setting it to ''random''.\n']);
        mosaic.noiseFlag = 'random';
    end

    if strcmp(mosaic.noiseFlag, 'none')
        noisyMRGCresponseInstances = [];
        return;
    end

    % Add noise via the class method
    noisyMRGCresponseInstances = mosaic.noisyResponseInstances( ...
        noiseFreeMRGCresponses, 'seed', noiseSeed);
end
