function [S_surround_LMS_norm, W_L_LMS_norm, W_M_LMS_norm, W_S_LMS_norm, stats] = normalize_LMS_surround_to_LM_v2(mLM, mLMS)
%NORMALIZE_LMS_SURROUND_TO_LM  Match LMS surround gains to an LM reference.
%
%   [S_surround_LMS_norm, W_L_LMS_norm, W_M_LMS_norm, W_S_LMS_norm, stats] = ...
%       normalize_LMS_surround_to_LM(mLM, mLMS)
%
%   Inputs
%       mLM   - LM-only mRGC mosaic (reference)
%       mLMS  - LMS mRGC mosaic (to be normalized)
%
%   Outputs
%       S_surround_LMS_norm - full LMS surround matrix after normalization
%                             (L, M, and S rows rescaled)
%       W_L_LMS_norm        - L-cone surround submatrix after scaling
%       W_M_LMS_norm        - M-cone surround submatrix after scaling
%       W_S_LMS_norm        - S-cone surround submatrix after scaling
%       stats               - struct with per-RGC L/M/S means, counts,
%                             and scale factors
%
%   Method
%       For each RGC, compute mean non-zero L- and M-cone surround weights
%       in both LM and LMS mosaics. Where both means are > 0, scale LMS
%       L and M surrounds so their means match the LM means. S surrounds
%       are rescaled using the same per-RGC factor as L. Zeros (no
%       connection) remain zeros.

%% 1. Grab LMS and LM surround matrices and size checks
S_LM  = mLM.rgcRFsurroundConeConnectivityMatrix;   % [nCones_LM  x nRGC]
S_LMS = mLMS.rgcRFsurroundConeConnectivityMatrix;  % [nCones_LMS x nRGC]

[~, nRGC_LM]   = size(S_LM);
[~, nRGC_LMS]  = size(S_LMS);

assert(nRGC_LM == nRGC_LMS, ...
    'LM and LMS mosaics must have the same number of RGCs.');

nRGC = nRGC_LM;

%% 2. Cone mosaics and cone-type indices
cm_LM  = mLM.inputConeMosaic;
cm_LMS = mLMS.inputConeMosaic;

lIdx_LM = cm_LM.lConeIndices;
mIdx_LM = cm_LM.mConeIndices;

lIdx_LMS = cm_LMS.lConeIndices;
mIdx_LMS = cm_LMS.mConeIndices;
sIdx_LMS = cm_LMS.sConeIndices;   

%% 3. Extract cone-type surround submatrices (L and M only)

% LM model (reference)
W_L_LM = full(S_LM(lIdx_LM, :));   % [nL_LM x nRGC]
W_M_LM = full(S_LM(mIdx_LM, :));   % [nM_LM x nRGC]

% LMS model (to be scaled)
W_L_LMS = full(S_LMS(lIdx_LMS, :)); % [nL_LMS x nRGC]
W_M_LMS = full(S_LMS(mIdx_LMS, :)); % [nM_LMS x nRGC]
W_S_LMS = full(S_LMS(sIdx_LMS, :)); % S-cone surround submatrix (LMS only)

%% 4. Per-RGC means (ignoring zeros) for LM and LMS

% ---- L-cone surrounds ----
Lnum_LM  = sum(W_L_LM  ~= 0, 1);   % # of non-zero L-cone weights (LM)
Lnum_LMS = sum(W_L_LMS ~= 0, 1);   % # of non-zero L-cone weights (LMS)

Lsum_LM  = sum(W_L_LM,  1);
Lsum_LMS = sum(W_L_LMS, 1);

Lmean_LM  = zeros(1, nRGC);
Lmean_LMS = zeros(1, nRGC);

idx_LM_nonzero  = (Lnum_LM  > 0);
idx_LMS_nonzero = (Lnum_LMS > 0);

Lmean_LM(idx_LM_nonzero)  = Lsum_LM(idx_LM_nonzero)  ./ Lnum_LM(idx_LM_nonzero);
Lmean_LMS(idx_LMS_nonzero) = Lsum_LMS(idx_LMS_nonzero) ./ Lnum_LMS(idx_LMS_nonzero);

% ---- M-cone surrounds ----
Mnum_LM  = sum(W_M_LM  ~= 0, 1);
Mnum_LMS = sum(W_M_LMS ~= 0, 1);

Msum_LM  = sum(W_M_LM,  1);
Msum_LMS = sum(W_M_LMS,  1);

Mmean_LM  = zeros(1, nRGC);
Mmean_LMS = zeros(1, nRGC);

idx_ML_nonzero  = (Mnum_LM  > 0);
idx_MLS_nonzero = (Mnum_LMS > 0);

Mmean_LM(idx_ML_nonzero)   = Msum_LM(idx_ML_nonzero)   ./ Mnum_LM(idx_ML_nonzero);
Mmean_LMS(idx_MLS_nonzero) = Msum_LMS(idx_MLS_nonzero) ./ Mnum_LMS(idx_MLS_nonzero);

% ---- S-cone surrounds (LMS only, pre-normalization) 
Snum_LMS = sum(W_S_LMS ~= 0, 1);      % # of non-zero S-cone weights (LMS)
Ssum_LMS = sum(W_S_LMS,  1);

Smean_LMS_in = zeros(1, nRGC);
idx_S_nonzero = (Snum_LMS > 0);
Smean_LMS_in(idx_S_nonzero) = Ssum_LMS(idx_S_nonzero) ./ Snum_LMS(idx_S_nonzero);

%% 5. Build per-RGC scale factors for L and M surrounds

% Only scale where BOTH LM and LMS means are > 0
scaleL = ones(1, nRGC);
maskL  = (Lmean_LM > 0) & (Lmean_LMS > 0);
scaleL(maskL) = Lmean_LM(maskL) ./ Lmean_LMS(maskL);

scaleM = ones(1, nRGC);
maskM  = (Mmean_LM > 0) & (Mmean_LMS > 0);
scaleM(maskM) = Mmean_LM(maskM) ./ Mmean_LMS(maskM);

%% 6. Apply scaling to LMS L- and M-cone surround weights

% Implicit expansion: each column j scaled by scaleL(j) or scaleM(j)
W_L_LMS_norm = W_L_LMS .* scaleL;   % L columns rescaled
W_M_LMS_norm = W_M_LMS .* scaleM;   % M columns rescaled

% Zeros remain zeros automatically.

%% 7. Reassemble normalized surround matrix for LMS (L and M rows replaced)

S_surround_LMS_norm = S_LMS;                % start from original LMS surround
S_surround_LMS_norm(lIdx_LMS, :) = W_L_LMS_norm;  % updated L rows
S_surround_LMS_norm(mIdx_LMS, :) = W_M_LMS_norm;  % updated M rows
% S-cone rows are unchanged.

%% 7b. Normalize S-cone surrounds using the same per-RGC scale 

scaleS       = scaleL;                       % reuse L-based per-RGC scale
W_S_LMS_norm = W_S_LMS .* scaleS;           % apply per-RGC scaling
S_surround_LMS_norm(sIdx_LMS, :) = W_S_LMS_norm;  % updated S rows

%% 8. Optional stats struct for inspection/debugging

stats.L.mean_LM       = Lmean_LM;
stats.L.mean_LMS_in   = Lmean_LMS;
stats.L.mean_LMS_out  = sum(W_L_LMS_norm, 1) ./ max(sum(W_L_LMS_norm ~= 0, 1), 1);
stats.L.num_LM        = Lnum_LM;
stats.L.num_LMS       = Lnum_LMS;
stats.L.scale         = scaleL;

stats.M.mean_LM       = Mmean_LM;
stats.M.mean_LMS_in   = Mmean_LMS;
stats.M.mean_LMS_out  = sum(W_M_LMS_norm, 1) ./ max(sum(W_M_LMS_norm ~= 0, 1), 1);
stats.M.num_LM        = Mnum_LM;
stats.M.num_LMS       = Mnum_LMS;
stats.M.scale         = scaleM;

% ---- S-cone surrounds (LMS only, post-normalization)
Snum_LMS_norm = sum(W_S_LMS_norm ~= 0, 1);
Ssum_LMS_norm = sum(W_S_LMS_norm, 1);

Smean_LMS_out = zeros(1, nRGC);
idx_S_norm_nonzero = (Snum_LMS_norm > 0);
Smean_LMS_out(idx_S_norm_nonzero) = Ssum_LMS_norm(idx_S_norm_nonzero) ./ Snum_LMS_norm(idx_S_norm_nonzero);

stats.S.mean_LMS_in   = Smean_LMS_in;
stats.S.mean_LMS_out  = Smean_LMS_out;
stats.S.num_LMS       = Snum_LMS;
stats.S.num_LMS_norm  = Snum_LMS_norm;
stats.S.scale         = scaleS;

%% 9. Package normalized surround info and save to disk
% This does NOT modify the mLMS object; it just saves the
% normalized surround connectivity and related stats in a helper struct.

LMSnorm = struct();
LMSnorm.S_surround = S_surround_LMS_norm;   % full normalized surround matrix
LMSnorm.W_L        = W_L_LMS_norm;          % normalized L-cone surround submatrix
LMSnorm.W_M        = W_M_LMS_norm;          % normalized M-cone surround submatrix
LMSnorm.W_S        = W_S_LMS_norm;          % normalized S-cone submatrix
LMSnorm.stats      = stats;

LMSnorm.description = [ ...
    'LMS surround normalized per RGC, separately for L and M, ' ...
    'and S-cone surrounds scaled with the same per-RGC factor as L.' ...
    ];

% Save to a helper file in the current directory
try
    save('mRGC_2deg_LMS_normBits.mat', 'LMSnorm', '-v7.3');
    fprintf('Saved LMSnorm surround weights to mRGC_2deg_LMS_normBits.mat\n');
catch ME
    warning(ME.identifier, 'Could not save LMSnorm struct: %s', ME.message);
end

end
