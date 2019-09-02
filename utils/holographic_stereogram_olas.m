%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Light field to hologram conversion using sliding-window-type holographic 
%   stereogram methods:
%       OLA-APAS - sliding window phase-added stereogram - overlap and add method (with depth)
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
%   Output: apas_ola    -   accurate phase-added stereogram (with depth)
%                           with sliding window (overlap and add) = 
%                           complex-valued wave field with phase and
%                           amplitude of size [Ny Nx C]
%
%   Example:    apas_ola = holographic_stereogram(lightField,lightFieldDepth);
%
%   Gordon Wetzstein
%   Stanford Computational Imaging Lab
%   gordon.wetzstein@stanford.edu
%   12/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [apas_ola] = holographic_stereogram_olas(lightField, lightFieldDepth, waveNum,...
    pixelPitch)

    if nargin < 3
        lambda = 532e-9; % 532 nm
        waveNum = [1 1 1] .* 2*pi/lambda;
    end
    
    if nargin < 4
        pixelPitch = 6.4e-6;
    end
    
    if length(waveNum) == 1
        waveNum = [1 1 1] * waveNum;
    end
    
    lambda = 2 * pi ./ waveNum;
            
    % hogel size is same as angular resolution
    hogelResolution = [size(lightField,1) size(lightField,2)];
    
    % resolution of hologram is same spatial resolution of light field
    hologramResolution = [size(lightField,3) size(lightField,4)];
    
    % number of color channels
    numColorChannels = 1;
    if ndims(lightField) > 4
        numColorChannels = size(lightField,5);
    end
    
    % initialize hologram with zeros, padded to avoid edges/for centering
    hogelRadius = floor(hogelResolution / 2);
    apas_ola = zeros([(hologramResolution + hogelRadius*2) numColorChannels]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute synthesis window

    % custom version of hann without zeros at ends
    function wndw = w_func(len)
        wndw = hann(len + 2);
        wndw = wndw(2:end-1);
    end

    win = w_func(hogelResolution(1)) * w_func(hogelResolution(1))';
    win = win / sum(win(:));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lfShape = size(lightField);

    % compute complex field
    complexLF = zeros(lfShape);
    for c = 1:numColorChannels

        % apply depth compensation
        compensatedDepth = zeros(size(lightFieldDepth));
        freqX = linspace(-1+1/hogelResolution(2), 1-1/hogelResolution(2), ...
            hogelResolution(2)) / (2 * pixelPitch);
        freqY = linspace(-1+1/hogelResolution(1), 1-1/hogelResolution(1), ...
            hogelResolution(1)) / (2 * pixelPitch);
        for ky = 1:hogelResolution(1)
            for kx = 1:hogelResolution(2)
                theta = asin(sqrt(freqX(kx)^2 + freqY(ky)^2) * lambda(c));
                compensatedDepth(ky, kx, :, :) = lightFieldDepth(ky, kx, :, :)...
                    * (1 - cos(theta));
            end
        end

        complexLF(:,:,:,:,c) = sqrt(lightField(:,:,:,:,c)) ...
            .* exp(1i * waveNum(c) .* compensatedDepth);
    end
    clear lightField lightFieldDepth compensatedDepth;

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

    % apply window
    complexLF = complexLF .* repmat(win, [1 1 lfShape(3:end)]);

    % overlap and add the hogels
    for ky = 1:hogelResolution(1)
        for kx = 1:hogelResolution(2)
            apas_ola(ky:ky+hologramResolution(1)-1, kx:kx+hologramResolution(2)-1, :) = ...
                apas_ola(ky:ky+hologramResolution(1)-1, kx:kx+hologramResolution(2)-1, :) ...
                + squeeze(complexLF(ky, kx, :, :, :));
        end
    end
    
    % crop back to light field size
    apas_ola = apas_ola(1+hogelRadius(1):end-hogelRadius(1),...
                        1+hogelRadius(2):end-hogelRadius(2), :);
end