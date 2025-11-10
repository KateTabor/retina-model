function scene = sceneFromMyXYZTiffv0(tiffPath, sceneFOVdeg)
%SCENEFROMMYXYZTIFF  Load our float32 XYZ (D65) TIFF and make an ISETBio scene.
%
%   scene = sceneFromMyXYZTiff('/path/to/initial_XYZ.tiff', 12);

    % 1) read the float32 XYZ tiff
    xyz = imread(tiffPath);   % for float32 3-ch TIFF this usually works
    xyz = double(xyz);        % HxWx3, D65-referred, 0..1-ish

    % (optional) check scale
    mx = max(xyz(:));
    if mx > 1.5
        warning('XYZ values look >1.5, consider dividing by 100 or checking writer.');
    end

    % 2) XYZ -> linear sRGB
    % needs isetcam/isetbio on path
    M_xyz2srgb = colorTransformMatrix('xyz2srgb');
    sz = size(xyz);
    rgb_lin = imageLinearTransform(reshape(xyz, [], 3), M_xyz2srgb);
    rgb_lin = reshape(rgb_lin, sz);

    % clip to 0..1
    rgb_lin = ieClip(rgb_lin, 0, 1);

    % 3) linear sRGB -> display sRGB (gamma)
    rgb_disp = lrgb2srgb(rgb_lin);  % ISETCam has this; if not, use ieLrgb2srgb

    % 4) make it uint8 so sceneFromFile is happy
    rgb_u8 = uint8(rgb_disp * 255);

    % 5) create an sRGB display
    d = displayCreate('LCD-Apple');

    % 6) finally: make the scene
    scene = sceneFromFile(rgb_u8, 'rgb', [], d);

    % 7) set field of view
    if nargin < 2
        sceneFOVdeg = 12;
    end
    scene = sceneSet(scene, 'fov', sceneFOVdeg);
end