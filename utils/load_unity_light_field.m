function [light_field, depth] = load_unity_light_field(datapath,...
    eyepieceFocalLength, frameNum, flipLFOutput, loadOnlyCentralView)

if nargin < 5
    loadOnlyCentralView = false;
end

% json calibration file name
json_fname = fullfile(datapath, 'cameras.json');

% read json file
json = jsondecode(fileread(json_fname));

% near clipping plane
zNear = json.NearClip;

% far clipping plane
zFar = json.FarClip;

% height and width of viewport plane
h = json.ViewportHeight;
w = json.ViewportWidth;

% get resolution and scaling factor for SLM units, forces all light fields
% to be same size if unitScale is set to imageWidth/ViewportWidth
slmPitch = 6.4e-6;
imageResolution = [json.PixelHeight, json.PixelWidth];
imageWidth = slmPitch * imageResolution(2);
if nargin >= 2 && ~isempty(eyepieceFocalLength)
    % scale imageWidth by magnification
    eyepieceVirtualImageDist = json.CameraDistance - eyepieceFocalLength;
    eyepieceHologramDist = 1/(1/eyepieceFocalLength + 1/eyepieceVirtualImageDist);
    magnification = eyepieceFocalLength / (eyepieceFocalLength - eyepieceHologramDist);
    imageWidth = imageWidth * magnification;
end
unitScale = imageWidth / json.ViewportWidth;

h = unitScale * h;
w = unitScale * w;
zNear = unitScale * zNear;
zFar = unitScale * zFar;

% get a grid for x and y coords in window coordinates
[xx_win, yy_win] = meshgrid(...
    linspace(0, imageResolution(2), imageResolution(2)),...
    linspace(imageResolution(1), 0, imageResolution(1)));
        
% calculate pixel positions given depth
xx_ndc = xx_win ./ imageResolution(2) - 1/2;
yy_ndc = yy_win ./ imageResolution(1) - 1/2;

if loadOnlyCentralView
    % specify the coordinates of the center view
    centerYView = floor(json.CameraRows/2) + 1;
    centerXView = floor(json.CameraColumns/2) + 1;
else
    % allocate memory for light field and depth
    light_field = zeros(json.CameraRows, json.CameraColumns,...
        imageResolution(1), imageResolution(2), 3);
    depth = zeros(json.CameraRows, json.CameraColumns,...
        imageResolution(1), imageResolution(2));
end

for camy=1:json.CameraRows
    % fprintf('%d', camy)
    % skip views if loading only central view
    if loadOnlyCentralView && camy ~= centerYView
        continue
    end

    for camx=1:json.CameraColumns
        % fprintf('.');
        % skip views if loading only central view
        if loadOnlyCentralView && camx ~= centerXView
            continue
        end

        % camera index, flip y coordinate
        camidx = (json.CameraRows-camy)*json.CameraColumns + camx;

        % camera position relative to central view
        campos_x = json.Cameras(camidx).parameters.localPosition.x;
        campos_y = json.Cameras(camidx).parameters.localPosition.y;

        campos_x = unitScale * campos_x;
        campos_y = unitScale * campos_y;

        % load depth map and light field view
        if nargin < 3 || isempty(frameNum)
            imageFilePath = fullfile(datapath,...
                [json.Cameras(camidx).key '_rgbd.png']);
        else
            imageFilePath = fullfile(datapath, json.Cameras(camidx).key,...
                sprintf('rgbd_%04d.png', frameNum));
        end
        [I, ~, D] = imread(imageFilePath);
        if ndims(I) == 2
            I = reshape([I I I], [size(I) 3]);
        end

        % convert to normalized double precision floating point values
        I = im2double(I);
        D = 1 ./ (im2double(D) .* (1 ./ zNear - 1 ./ zFar) + 1 ./ zFar);

        % get/reset zero disparity plane
        zero_disp_plane = json.CameraDistance;
        zero_disp_plane = unitScale * zero_disp_plane;

        % target position on SLM/viewport/zero_disparity_plane for each pixel
        xx_slm = xx_ndc * w;
        yy_slm = yy_ndc * h;

        % account for camera position's depth-dependent shift
        x_offset = (zero_disp_plane - D) / zero_disp_plane * campos_x;
        y_offset = (zero_disp_plane - D) / zero_disp_plane * campos_y;

        % point cloud relative to central camera position
        xx_metric = xx_ndc * w .* D / zero_disp_plane + x_offset;
        yy_metric = yy_ndc * h .* D / zero_disp_plane + y_offset;

        % use focal length to convert depth to be relative to hologram
        % plane (which is assumed to be the zero disparity plane)
        if nargin >= 2 && ~isempty(eyepieceFocalLength)
            virtualImageDist = D - eyepieceFocalLength;
            imageDist = 1 ./ (1/eyepieceFocalLength + 1./virtualImageDist);
            imageMag = eyepieceFocalLength ./ (eyepieceFocalLength - imageDist);
            virtualZeroDisp = zero_disp_plane - eyepieceFocalLength;

            zero_disp_plane = 1 ./ (1/eyepieceFocalLength + 1./virtualZeroDisp);
            zeroDispMag = eyepieceFocalLength ./ (eyepieceFocalLength - zero_disp_plane);

            xx_metric = xx_metric ./ imageMag;
            yy_metric = yy_metric ./ imageMag;
            xx_slm = xx_slm ./ zeroDispMag;
            yy_slm = yy_slm ./ zeroDispMag;
            D = imageDist;
        end

        % positions relative to corresponding SLM pixel
        xx_dist = xx_slm - xx_metric;
        yy_dist = yy_slm - yy_metric;
        zz_dist = zero_disp_plane - D;

        % distance from pixel to corresponding SLM pixel
        abs_dist = sqrt(xx_dist.^2 + yy_dist.^2 + zz_dist.^2);
        % sign for which side of slm
        metric_dist = abs_dist .* sign(zz_dist);

        if loadOnlyCentralView
            light_field = I;
            depth = metric_dist;
        else
            light_field(camy, camx, :, :, :) = I;
            depth(camy, camx, :, :) = metric_dist;
        end

    end
end

if nargin >= 4 && flipLFOutput && ~loadOnlyCentralView
    light_field = flip(light_field, 1);
    light_field = flip(light_field, 2);
    depth = flip(depth, 1);
    depth = flip(depth, 2);
    depth = -depth;
end


end