function theMRGCmosaic_S = build2degONmRGCMosaic_enableSconeSurrounds(varargin)
% build2degONmRGCMosaic_enableSconeSurrounds
% Enable S-cone input to the surrounds of ALL mRGCs and save the result.
%
% Syntax:
%   theMRGCmosaic_S = build2degONmRGCMosaic_enableSconeSurrounds
%   theMRGCmosaic_S = build2degONmRGCMosaic_enableSconeSurrounds('opticsSubjectName','PLOSpaperDefaultSubject', ...)
%
% Name-Value pairs (match your v2 defaults):
%   'coneMosaicSpecies'         - 'human' (default) or 'macaque'
%   'opticsSubjectName'         - default: 'PLOSpaperDefaultSubject'
%   'rgcMosaicName'             - default: 'PLOSpaperTemporal2DegsMosaic'
%   'targetVisualSTFdescriptor' - default: 'default'
%   'indicesOfSconeSurroundPoolingEnabledRGCs' - [] (default = ALL RGCs)
%   'exportDir'                 - output folder (default: 'synthetic_retina_exports')
%   'exportMatFile'             - true/false (default: true)
%   'closeOpenFigures'          - true/false (default: true)
%
% Returns:
%   theMRGCmosaic_S  - the modified mRGC mosaic with S-cone surrounds enabled
%
% Notes:
% - Uses the same underlying cone/mRGC positions as your v2 when you pass the same 4 spec args.
% - Centers remain unchanged; only surround pooling is modified.

% --- keep caller args safe
userArgs = varargin;

% --- init ISETBio/ISETCam safely (avoid varargin being cleared)
isInit = false;
if exist('ieSessionGet','file') == 2
    try isInit = ieSessionGet('initialized'); catch, isInit = false; end
end
if ~isInit && exist('ieInit','file') == 2
    try evalin('base','ieInit;'); catch, ieInit; end
end

% --- parse inputs
p = inputParser;
p.addParameter('coneMosaicSpecies','human');
p.addParameter('opticsSubjectName','PLOSpaperDefaultSubject');
p.addParameter('rgcMosaicName','PLOSpaperTemporal2DegsMosaic');
p.addParameter('targetVisualSTFdescriptor','default');
p.addParameter('indicesOfSconeSurroundPoolingEnabledRGCs',[]);
p.addParameter('exportDir','synthetic_retina_exports');
p.addParameter('exportMatFile',true);
p.addParameter('closeOpenFigures',true);
p.parse(userArgs{:});
opts = p.Results;

if opts.closeOpenFigures, close all; end

% --- load the SAME prebaked mosaic spec as v2 (optics not required here)
[theMRGCmosaic_S, ~, ~] = mRGCMosaic.loadPrebakedMosaic( ...
    opts.coneMosaicSpecies, ...
    opts.opticsSubjectName, ...
    opts.rgcMosaicName, ...
    opts.targetVisualSTFdescriptor, ...
    'computeTheMosaicOptics', false);

% --- choose target RGCs (default = ALL)
if isempty(opts.indicesOfSconeSurroundPoolingEnabledRGCs)
    targetRGCs = 1:theMRGCmosaic_S.rgcsNum;
else
    targetRGCs = opts.indicesOfSconeSurroundPoolingEnabledRGCs(:)';
end

fprintf('Enabling S-cone surround pooling for %d / %d mRGCs ...\n', ...
    numel(targetRGCs), theMRGCmosaic_S.rgcsNum);

% --- apply modification (no mass visualization to avoid huge output)
theMRGCmosaic_S.enableSconeSurroundPoolingInSelectCells( ...
    targetRGCs, ...
    'visualizeBeforeAndAfterConePoolingMaps', false, ...
    'exportVisualizationPDFdirectory', '' );

% --- summarize
fprintf('Done. Centers unchanged; surrounds now pool indiscriminately from L/M/S.\n');

% % --- optional export for future reuse
% if opts.exportMatFile
%     if ~exist(opts.exportDir,'dir'), mkdir(opts.exportDir); end
% 
%     % Build a simple, informative filename from spec + mosaic geometry
%     eccStr  = sprintf('Ecc_%0.1f_%0.1f', theMRGCmosaic_S.eccentricityDegs(1), theMRGCmosaic_S.eccentricityDegs(2));
%     sizeStr = sprintf('Size_%0.1fx%0.1f', theMRGCmosaic_S.sizeDegs(1), theMRGCmosaic_S.sizeDegs(2));
%     tag     = sprintf('%s_%s_%s_%s_%s_%s_SconeSurrounds', ...
%         opts.coneMosaicSpecies, opts.opticsSubjectName, opts.rgcMosaicName, ...
%         opts.targetVisualSTFdescriptor, eccStr, sizeStr);
%     tag = regexprep(tag,'[^A-Za-z0-9_\.]+',''); % sanitize
% 
%     outPath = fullfile(opts.exportDir, ['MRGCMosaic_' tag '.mat']);
%     theMRGCmosaic = theMRGCmosaic_S; %#ok<NASGU>  % save under a familiar name
%     metadata = struct('generatedOn', datestr(now), 'spec', opts, ...
%                       'note','S-cone surround pooling enabled for the listed RGC indices'); %#ok<NASGU>
%     save(outPath, 'theMRGCmosaic', 'metadata','-v7.3');
%     fprintf('Saved modified mosaic to:\n  %s\n', outPath);
% end
% end
