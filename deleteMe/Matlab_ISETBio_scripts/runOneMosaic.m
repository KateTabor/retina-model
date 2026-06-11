function [respMean, coneRespSingle, rgcTimeSeries] = runOneMosaic(mosaic, oi, dt, nFrames, timeAxis, surroundOverride)
% RUNONEMOSAIC  Deterministic mRGC run + optional cone/mRGC outputs.
%
%   respMean = runOneMosaic(mosaic, oi, dt, nFrames, timeAxis)
%       -> mean mRGC response over time, [1 x Nrgc]
%
%   [respMean, coneRespSingle] = runOneMosaic(...)
%       -> also return single-frame cone response to this OI
%
%   [respMean, coneRespSingle, rgcTimeSeries] = runOneMosaic(...)
%       -> also return full time series of mRGC responses [T x Nrgc]
%
%   NEW (Nov 2025):
%     If a 6th argument `surroundOverride` is provided and non-empty, this
%     function calls compute_mRGC_withCustomSurround(), which uses that
%     matrix as the surround connectivity instead of the mosaic's internal
%     rgcRFsurroundConeConnectivityMatrix.

    % ----- NEW: handle optional surroundOverride argument -----
    if nargin < 6                                      % NEW
        surroundOverride = [];                         % NEW
    end                                                % NEW
    useCustomSurround = ~isempty(surroundOverride);    % NEW

    % Default outputs
    coneRespSingle = [];
    rgcTimeSeries  = [];

    % --- Cone mosaic feeding this mRGC mosaic
    cm = mosaic.inputConeMosaic;

    % Cones: deterministic integration
    cm.noiseFlag = 'frozen';
    if isprop(cm, 'integrationTime')
        cm.integrationTime = dt;
    end

    % mRGC noise: also deterministic
    if isprop(mosaic, 'noiseFlag')
        mosaic.noiseFlag = 'frozen';
    end

    % --- Single-frame cone response to the optical image
    coneResp = cm.compute(oi);
    coneResp = double(squeeze(coneResp));

    % --- Reshape to [1 x T x Ncones] and replicate over time
    if isvector(coneResp)
        nCones = numel(coneResp);
        cone3D = reshape(coneResp, [1 1 nCones]);
    else
        sz = size(coneResp);
        if numel(sz)==2 && all(sz > 1)
            if sz(1) > sz(2)
                tmp = coneResp.';
            else
                tmp = coneResp;
            end
            cone3D = reshape(tmp, [1 size(tmp,1) size(tmp,2)]);
        else
            if numel(sz) < 3
                coneResp(:,:,2) = coneResp;
                sz = size(coneResp);
            end
            tmp    = reshape(coneResp, [], sz(3));
            cone3D = reshape(tmp.', [1 sz(3) numel(tmp)/sz(3)]);
        end
    end

    if size(cone3D,2) == 1
        cone3D = repmat(cone3D, [1 nFrames 1]);
    else
        if size(cone3D,2) < nFrames
            cone3D = repmat(cone3D, [1 ceil(nFrames/size(cone3D,2)) 1]);
        end
        cone3D = cone3D(:,1:nFrames,:);
    end

    % --- Run mRGC mosaic (deterministic seed)
    if useCustomSurround
        % NEW: call our wrapper that overrides surround weights
        [rgcResp, rgcRespNoisy, ~] = compute_mRGC_withCustomSurround( ...    % NEW
            mosaic, cone3D, timeAxis, surroundOverride, 'seed', 1);          % NEW
    else
        % Original behavior: call the built-in compute method
        [rgcResp, rgcRespNoisy] = mosaic.compute(cone3D, timeAxis, 'seed', 1);  % UNCHANGED
    end

    % Prefer noisy instances if present (with frozen noise it's deterministic)
    if ~isempty(rgcRespNoisy)
        X = squeeze(rgcRespNoisy);   % [T x Nrgc]
    else
        X = squeeze(rgcResp);        % [T x Nrgc]
    end
    if isvector(X)
        X = X(:)';                   % [1 x Nrgc]
    end

    % Mean over time
    respMean = mean(X,1);            % [1 x Nrgc]

    % Optional extra outputs
    if nargout > 1
        coneRespSingle = coneResp;   % single-frame cone response
    end
    if nargout > 2
        rgcTimeSeries = X;           % [T x Nrgc] time series
    end
end




% function [respMean, coneRespSingle, rgcTimeSeries] = runOneMosaic(mosaic, oi, dt, nFrames, timeAxis)
% % RUNONEMOSAIC  Deterministic mRGC run + optional cone/mRGC outputs.
% %
% %   respMean = runOneMosaic(mosaic, oi, dt, nFrames, timeAxis)
% %       -> mean mRGC response over time, [1 x Nrgc]
% %
% %   [respMean, coneRespSingle] = runOneMosaic(...)
% %       -> also return single-frame cone response to this OI
% %
% %   [respMean, coneRespSingle, rgcTimeSeries] = runOneMosaic(...)
% %       -> also return full time series of mRGC responses [T x Nrgc]
% 
%     % Default outputs
%     coneRespSingle = [];
%     rgcTimeSeries  = [];
% 
%     % --- Cone mosaic feeding this mRGC mosaic
%     cm = mosaic.inputConeMosaic;
% 
%     % Cones: deterministic integration
%     cm.noiseFlag = 'frozen';
%     if isprop(cm, 'integrationTime')
%         cm.integrationTime = dt;
%     end
% 
%     % mRGC noise: also deterministic
%     if isprop(mosaic, 'noiseFlag')
%         mosaic.noiseFlag = 'frozen';
%     end
% 
%     % --- Single-frame cone response to the optical image
%     coneResp = cm.compute(oi);
%     coneResp = double(squeeze(coneResp));
% 
%     % --- Reshape to [1 x T x Ncones] and replicate over time
%     if isvector(coneResp)
%         nCones = numel(coneResp);
%         cone3D = reshape(coneResp, [1 1 nCones]);
%     else
%         sz = size(coneResp);
%         if numel(sz)==2 && all(sz > 1)
%             if sz(1) > sz(2)
%                 tmp = coneResp.';
%             else
%                 tmp = coneResp;
%             end
%             cone3D = reshape(tmp, [1 size(tmp,1) size(tmp,2)]);
%         else
%             if numel(sz) < 3
%                 coneResp(:,:,2) = coneResp;
%                 sz = size(coneResp);
%             end
%             tmp    = reshape(coneResp, [], sz(3));
%             cone3D = reshape(tmp.', [1 sz(3) numel(tmp)/sz(3)]);
%         end
%     end
% 
%     if size(cone3D,2) == 1
%         cone3D = repmat(cone3D, [1 nFrames 1]);
%     else
%         if size(cone3D,2) < nFrames
%             cone3D = repmat(cone3D, [1 ceil(nFrames/size(cone3D,2)) 1]);
%         end
%         cone3D = cone3D(:,1:nFrames,:);
%     end
% 
%     % --- Run mRGC mosaic (deterministic seed)
%     [rgcResp, rgcRespNoisy] = mosaic.compute(cone3D, timeAxis, 'seed', 1);
% 
%     % Prefer noisy instances if present (with frozen noise it's deterministic)
%     if ~isempty(rgcRespNoisy)
%         X = squeeze(rgcRespNoisy);   % [T x Nrgc]
%     else
%         X = squeeze(rgcResp);        % [T x Nrgc]
%     end
%     if isvector(X)
%         X = X(:)';                   % [1 x Nrgc]
%     end
% 
%     % Mean over time
%     respMean = mean(X,1);            % [1 x Nrgc]
% 
%     % Optional extra outputs
%     if nargout > 1
%         coneRespSingle = coneResp;   % single-frame cone response
%     end
%     if nargout > 2
%         rgcTimeSeries = X;           % [T x Nrgc] time series
%     end
% end
