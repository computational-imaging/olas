%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Light field to hologram conversion using the holographic stereogram
%   method or, if depth for each light ray is available, the phase-added 
%   stereogram method.
%
%   Note:   if no depth information is specified, the algorithm will compute
%           the holographic stereogram (HS) and with depth data it will 
%           compute the phase-added stereogram (PAS).
%
%   Input:  lightField -    light field with size [My Mx Ny Nx C]       
%                               My, Mx are number of angular samples
%                               Ny, Nx are number of spatial samples
%                               C are number of color channels
%
%           lightFieldDepth -   optional light field with metric depth values per ray with size [My Mx Ny Nx]       
%                               My, Mx are number of angular samples
%                               Ny, Nx are number of spatial samples
%
%           k               -   wavenumbers, i.e. 2*pi/lambda for each
%                               color channel (only need to specify when
%                               using with depth)
%
%   Output: hs -    holographic stereogram = complex-valued wave field with phase and amplitude of size [Ny Nx C]
%           pas -   phase-added stereogram = complex-valued wave field with phase and amplitude of size [Ny Nx C]
%
%   Example:    hs = holographic_stereogram(lightField);
%               [hs,pas] = holographic_stereogram(lightField,lightFieldDepth,k);
%
%   Gordon Wetzstein
%   Stanford Computational Imaging Lab
%   gordon.wetzstein@stanford.edu
%   12/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hs] = holographic_stereogram_hs(lightField, ~, k)

    if nargin < 3
        %disp('WARNING: need to specify wavenumber! Using default values.');
        lambda = 532e-9; % 532 nm
        k = [1 1 1] .* 2*pi/lambda;
    end

    if length(k) == 1
        k = [1 1 1] * k;
    end
                
    % hogel size is same as angular resolution
    hogelResolution = [size(lightField,1) size(lightField,2)];

    % target resolution of hologram is same spatial resolution of light field
    hologramCroppedResolution = [size(lightField,3) size(lightField,4)];
    
    % pad so that we fit full hogels of the next size up
    numHogels = ceil(hologramCroppedResolution ./ hogelResolution);
    hologramCropAmt = numHogels .* hogelResolution - hologramCroppedResolution;
    % make the crop amount a multiple of 2 for centered light field padding
    cropRadius = [ceil(hologramCropAmt / 2); floor(hologramCropAmt / 2)];
    hologramResolution = numHogels .* hogelResolution;

    % pad light field to match hologram
    lightField = padarray(lightField, [0 0 cropRadius(1,:)], 'pre');
    lightField = padarray(lightField, [0 0 cropRadius(2,:)], 'post');
    
    % number of color channels
    numColorChannels = 1;
    if ndims(lightField) > 4
        numColorChannels = size(lightField,5);
    end
    
    % center pixel within each hogel
    hogelCenter = floor(hogelResolution/2)+1;

    % number of hogels over SLM
    numTiledHogels = floor(hologramResolution./hogelResolution);
           
    % initialize hologram with zeros
    hs = zeros([hologramResolution numColorChannels]); 
    
    hc = hogelCenter; hr = hogelResolution;
    complexLF = sqrt(lightField(:,:,hc(1):hr(1):end,hc(2):hr(2):end,:));
    lfShape = size(complexLF);

    % fft the first dimension
    complexLF = reshape(complexLF, lfShape(1), []);
    complexLF = fftshift(fft(ifftshift(complexLF, 1), [], 1), 1) / lfShape(1);
    complexLF = reshape(complexLF, lfShape);
    % fft the second dimension
    complexLF = permute(complexLF, [2 1 3 4 5]);
    complexLF = reshape(complexLF, lfShape(2), []);
    complexLF = fftshift(fft(ifftshift(complexLF, 1), [], 1), 1) / lfShape(2);
    complexLF = reshape(complexLF, [lfShape(2) lfShape(1) lfShape(3:end)]);
    complexLF = permute(complexLF, [2 1 3 4 5]);
    
    % tile the hogels
    for ky = 1:hogelResolution(1)
        for kx = 1:hogelResolution(2)
            hs(ky:hogelResolution(1):end, kx:hogelResolution(2):end, : ) =...
                squeeze(complexLF(ky, kx, :, :, :));
        end
    end

    % crop back to original light field size
    hs = hs(1+cropRadius(1,1):end-cropRadius(2,1),...
            1+cropRadius(1,2):end-cropRadius(2,2),:);

    
end