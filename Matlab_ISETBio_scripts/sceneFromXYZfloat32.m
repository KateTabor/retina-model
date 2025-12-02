function scene = sceneFromXYZfloat32(tiffPath, sceneFOVdeg)
% sceneFromXYZfloat32  Load float32 XYZ (D65) TIFF and make an ISETBio scene.
%
%   scene = sceneFromXYZfloat32('/path/to/image_XYZ.tiff', 2);

    % --- 1) Read XYZ image ---
    xyz = double(imread(tiffPath));   % H x W x 3, float32 -> double

    if size(xyz,3) ~= 3
        error('Expected HxWx3 XYZ image. Got size: %s', mat2str(size(xyz)));
    end

    if max(xyz(:)) > 2
        warning('XYZ seems >2.0; confirm writer scaling (expected ~0..1).');
    end

    % --- 2) XYZ -> linear sRGB ---
    if exist('colorTransformMatrix','file') == 2
        M_xyz2srgb = colorTransformMatrix('xyz2srgb');
    else
        % Fallback matrix (standard XYZ->sRGB)
        M_xyz2srgb = [ 3.2406 -1.5372 -0.4986; ...
                      -0.9689  1.8758  0.0415; ...
                       0.0557 -0.2040  1.0570 ];
    end

    sz      = size(xyz);
    rgb_lin = imageLinearTransform(reshape(xyz, [], 3), M_xyz2srgb);
    rgb_lin = reshape(rgb_lin, sz);
    rgb_lin = ieClip(rgb_lin, 0, 1);

    % --- 3) linear sRGB -> display sRGB (gamma-encoded), still float ---
    if exist('lrgb2srgb','file') == 2
        rgb_disp = lrgb2srgb(rgb_lin);
    else
        % Fallback sRGB EOTF
        a = 0.055;
        rgb_disp = max(rgb_lin, 0);
        small = rgb_disp <= 0.0031308;
        rgb_disp(small)  = 12.92 * rgb_disp(small);
        rgb_disp(~small) = (1+a) .* rgb_disp(~small).^(1/2.4) - a;
    end
    rgb_disp = ieClip(rgb_disp, 0, 1);

    % --- 4) Build spectral scene using a calibrated display ---
    d     = displayCreate('LCD-Apple');
    scene = sceneFromFile(rgb_disp, 'rgb', [], d);

    % --- 5) Set field of view ---
    if nargin < 2 || isempty(sceneFOVdeg)
        sceneFOVdeg = 2;
    end
    scene = sceneSet(scene, 'fov', sceneFOVdeg);
end
