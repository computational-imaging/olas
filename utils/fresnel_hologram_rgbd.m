%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   load a light field + depth, keep only the central view, compte a
%   fresnel hologram from it by propagating every point separately
%       ATTENTION: very slow!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hologram] = fresnel_hologram_rgbd(rgb, depth, wavenums, ...
	phaseCompensationForPropagation, pixelPitch)

    if nargin < 3
        lambdas = 532e-9 * [1 1 1]; % 532 nm
    else
		lambdas = 2 * pi ./ wavenums;
	end
 
	% phase compensation on or off
	if nargin < 4
		phaseCompensationForPropagation = true;
	end
    
    if nargin < 5
        pixelPitch = 6.4e-6;
    end
	% SLM pixel size in m
	dx = pixelPitch; dy = pixelPitch;

	% list all color channels you want to process here
	channels = 1:size(rgb, 3);

	% downsample to reduce time
	imageResolution = [1080 1920]./2;

	% resize and convert to amplitude
	rgb = sqrt(imresize(rgb, imageResolution));
	depth = imresize(depth, imageResolution);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% compute Fresnel hologram

	% initialize hologram
	hologram = zeros([size(rgb,1) size(rgb,2) numel(channels)]);

	for kc=1:numel(channels)
	    
	    % current color channel
	    c = channels(kc);
	    
	    for ky=1:size(rgb,1)

	        disp(['Processing hologram, channel: ' num2str(kc) ' | ' num2str(numel(channels)) ', line ' num2str(ky) ' | ' num2str(size(rgb,1))]);

	        for kx=1:size(rgb,2)

	            % get current depth value
	            currentDepth = depth(ky,kx);

	            % compute propagation kernel for this depth
	            kernel = propagation_kernel(currentDepth, lambdas(kc), [dy dx], phaseCompensationForPropagation);

	            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	            % get indices to splat the kernel to
	            
	            cutTop      = max(1 - (ky-floor(size(kernel,1)/2)),0);
	            cutBottom   = max((ky+floor(size(kernel,1)/2)) - size(rgb,1), 0);
	            
	            cutLeft     = max(1 - (kx-floor(size(kernel,2)/2)),0);
	            cutRight    = max((kx+floor(size(kernel,2)/2)) - size(rgb,2), 0);
	            
	            % these are the coordinates for the kernel            
	            kernelIdxY = 1:size(kernel,1); kernelIdxY = kernelIdxY(1+cutTop:end-cutBottom);                        
	            kernelIdxX = 1:size(kernel,2); kernelIdxX = kernelIdxX(1+cutLeft:end-cutRight);
	            % get x and y indices
	            [idxXkernel, idxYkernel] = meshgrid(kernelIdxX,kernelIdxY);
	            
	            % these are the coordinates for the image
	            startY  = max(1,ky-floor(size(kernel,1)/2));
	            endY    = min(size(rgb,1),ky+floor(size(kernel,1)/2));            
	            startX  = max(1,kx-floor(size(kernel,2)/2));
	            endX    = min(size(rgb,2),kx+floor(size(kernel,2)/2));
	            % get x and y indices
	            [idxX, idxY] = meshgrid(startX:endX,startY:endY);                        
	            
	            % convert to matrix indices
	            linearIndImg    = sub2ind(size(depth), idxY(:), idxX(:));
	            linearIndKernel = sub2ind(size(kernel), idxYkernel(:), idxXkernel(:));
	            
	            
	            % temp array
	            hologramTmp = zeros([size(rgb,1) size(rgb,2)]);
	            hologramTmp(linearIndImg) = rgb(ky,kx,c) .* kernel(linearIndKernel);
	            
	            % splat into hologram
	            hologram(:,:,kc) = hologram(:,:,kc) + hologramTmp;
	            
	        end
	        
	        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	        % plot intermediate results

	        subplot(1,2,1);
	        imagesc(abs(hologram));

	        subplot(1,2,2);
	        imagesc(angle(hologram));
	        drawnow;
	    end
	end

end