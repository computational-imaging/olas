%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Light field to hologram conversion using tiled-window-type holographic 
%   stereogram method:
%       APAS    - sliding window phase-added stereogram - accurate methodm
%                  i.e. only keep window center (with depth)
%
%   Input:  lightField -    light field with size [My Mx Ny Nx C]       
%                               My, Mx are number of angular samples
%                               Ny, Nx are number of spatial samples
%                               C are number of color channels
%
%           lightFieldDepth -   optional light field with metric depth
%                               values per ray with size [My Mx Ny Nx]       
%                               My, Mx are number of angular samples
%                               Ny, Nx are number of spatial samples
%
%           waveNum         -   wavenumbers, i.e. 2*pi/lambda for each
%                               color channel (only need to specify when
%                               using with depth)
%
%   Output: apas        -   accurate phase-added stereogram (with depth)
%                           with sliding window, but only keep center
%                           value = complex-valued wave field with phase
%                           and amplitude of size [Ny Nx C]
%           
%   Example:    apas = holographic_stereogram(lightField,lightFieldDepth);
%
%   Gordon Wetzstein
%   Stanford Computational Imaging Lab
%   gordon.wetzstein@stanford.edu
%   12/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [apas] = holographic_stereogram_apas(lightField, lightFieldDepth, waveNum)

    if nargin < 3
        lambda = 532e-9; % 532 nm
        waveNum = [1 1 1] .* 2*pi/lambda;
    end
    
    if length(waveNum) == 1
        waveNum = [1 1 1] * waveNum;
    end
            
    % hogel size is set to 1/3 of the angular resolution
    outerHogelResolution = [size(lightField,1) size(lightField,2)];
    hogelResolution = ceil(outerHogelResolution / 3);
    hogelCropPx = (outerHogelResolution - hogelResolution) / 2;
    
    % resolution of hologram is same spatial resolution of light field
    hologramResolution = [size(lightField,3) size(lightField,4)];

    % number of hogels and indexing offsets
    hogelRadius = floor(hogelResolution / 2);
    numHogels = floor(hologramResolution ./ hogelResolution);
    hogelOffset = floor((hologramResolution - numHogels .* hogelResolution) / 2);
    
    % number of color channels
    numColorChannels = 1;
    if ndims(lightField) > 4
        numColorChannels = size(lightField,5);
    end
    
    % initialize hologram with zeros
    apas = zeros([hologramResolution numColorChannels]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 0:numHogels(1)-1
        % indices are corners of the hogels
        ky = (i * hogelResolution(1)) + hogelOffset(1) + 1;
        for j = 0:numHogels(2)-1
            kx = (j * hogelResolution(2)) + hogelOffset(2) + 1;

            for c=1:numColorChannels

                % get angular samples of light field and depth at idx, with
                % offset by hogelRadius to get the center pixel spatially
                % for correct angular neighborhood
                lfPatch = squeeze(lightField(:,:,ky+hogelRadius(1),kx+hogelRadius(2),c));
                depPatch = squeeze(lightFieldDepth(:,:,ky+hogelRadius(1),kx+hogelRadius(2)));

                % compute hogel via fft       
                hogel = fftshift(fft2(ifftshift( sqrt(lfPatch) ...
                    .* exp( 1i.*waveNum(c).*depPatch) ))) ./ prod(hogelResolution);  % complex-valued

                % crop and store center of hogel
                apas(ky:ky+hogelResolution(1)-1, kx:kx+hogelResolution(2)-1,c) = ...
                    hogel(1+hogelCropPx(1):end-hogelCropPx(1),...
                          1+hogelCropPx(2):end-hogelCropPx(2));
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
    end

end