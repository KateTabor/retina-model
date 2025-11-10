function oi = makeOptics_humanWVF(varargin)
% makeOptics_humanWVF
% Create a reusable human wavefront-based OI (scene-agnostic), with
% optional custom wavelength sampling, provenance metadata, and
% optional saving.
%
% Examples:
%   oi = makeOptics_humanWVF;                              % default pupil, default WVF waves
%   oi = makeOptics_humanWVF('pupilDiameterMM', 2.5);      % set pupil
%   oi = makeOptics_humanWVF('wave', 400:5:700);           % custom wavelength grid
%   oi = makeOptics_humanWVF('saveToFile', true);          % save to auto-named file (-v7)
%   oi = makeOptics_humanWVF('saveToFile', true, ...
%                            'outFile','my_optics.mat','saveFormat','-v7.3');

% --- parse inputs (simple)
p = inputParser;
p.addParameter('pupilDiameterMM', 3.0, @isscalar);
p.addParameter('wave', [], @(x)(isnumeric(x)&&isvector(x)));
p.addParameter('saveToFile', false, @islogical);
p.addParameter('outFile','', @(s)(ischar(s)||isstring(s)));
p.addParameter('saveFormat','-v7', @(s) any(strcmp(s,{'-v7','-v7.3'})));
p.parse(varargin{:});
opts = p.Results;

% --- Init ISETBio/ISETCam (your exact block)
isInit = false;
if exist('ieSessionGet','file') == 2
    try, isInit = ieSessionGet('initialized'); catch, isInit = false; end
end
if ~isInit && exist('ieInit','file') == 2
    try, evalin('base','ieInit;'); catch, ieInit; end
end

% --- Start from human WVF optics OI, then modify WVF as needed
oi  = oiCreate('wvf human');           % <-- correct way to get human/Thibos WVF
wvf = oiGet(oi,'optics wvf');          % pull the WVF out

% Optional custom wavelength sampling (set BEFORE compute)
if ~isempty(opts.wave)
    wvf = wvfSet(wvf,'calc wavelengths', opts.wave(:)');
end

% Set pupil and compute WVF
wvf = wvfSet(wvf,'calc pupil diameter', opts.pupilDiameterMM, 'mm');
wvf = wvfCompute(wvf);

% Rebuild OI from WVF
oi = wvf2oi(wvf);

% --- Name + minimal provenance metadata
w = oiGet(oi,'wave');
if ~isempty(opts.wave)
    step = round(median(diff(w)));
    waveTag = sprintf('_W%dto%d_%dnm', w(1), w(end), step);
else
    waveTag = '';
end
oiName = sprintf('HumanWVF_%0.1fmmPupil%s', opts.pupilDiameterMM, waveTag);
oi = oiSet(oi,'name',oiName);

meta = struct();
meta.createdOn       = datestr(now);
meta.pupilDiameterMM = opts.pupilDiameterMM;
meta.wave_nm         = w(:);
meta.notes           = 'Scene-agnostic human WVF optics; FOV comes from scene.';
oi = oiSet(oi,'metadata',meta);

fprintf('Created OI: %s | pupil=%.2f mm | wave=[%d..%d] nm (N=%d)\n', ...
    oiName, opts.pupilDiameterMM, w(1), w(end), numel(w));

% --- Optional save
if opts.saveToFile
    outFile = char(opts.outFile);
    if isempty(strtrim(outFile))
        outFile = sprintf('optics_humanWVF_%0.1fmm%s.mat', opts.pupilDiameterMM, waveTag);
        outFile = strrep(outFile,'..','.'); % tidy
    end
    save(outFile,'oi',opts.saveFormat);
    fprintf('Saved optics to: %s (%s)\n', outFile, opts.saveFormat);
end
end
