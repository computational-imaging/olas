%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the propagation kernel for free-space propagation
%
%   input:  distance        - distance from point to hologram plane
%           lambda          - wavelength
%           slmPixelSize    - feature size of hologram
%           bPhaseCompensation - subtract distance normal to hologram plane
%                                   (i.e. "distance") from distance to each pixel

function h = propagation_kernel(distance, lambda, slmPixelSize, bPhaseCompensation)

    if nargin<4
        bPhaseCompensation = true;
    end
    
    % wave number
    k = 2*pi/lambda;
    
    % check if we're propagating back to SLM plane
    bBackpropagation = true;
    if distance>=0
        bBackpropagation = false;
    end
    distance = abs(distance);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % compute size of subhologram via grating equation 
    % d * (sin(theta_in) + sin(theta_out)) = lambda
    % and then use similar triangles with propagation distance       
    
    % max diffraction angle 
    maxDiffractionAngleRad = asin(lambda ./ (2*slmPixelSize));
    
    % size of subhologram in m
    subhologramSize = 2 * distance * tan(maxDiffractionAngleRad);
        
    % resolution of subhologram in number of slm pixels
    subhologramResolution = round((subhologramSize ./ slmPixelSize) );
    
    % this should be even
    %subhologramResolution = floor(subhologramResolution./2) .* 2 ;
    subhologramResolution = floor(subhologramResolution./2) .* 2 + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % compute convolution kernel
    
    % get coordinates of pixel centers in kernel
    startY = -slmPixelSize(1)*(subhologramResolution(1)-1)/2;
    startX = -slmPixelSize(2)*(subhologramResolution(2)-1)/2;
    [XKERNEL, YKERNEL] = meshgrid(   startX + ((0:subhologramResolution(2)-1)*slmPixelSize(2)), ...
                                    startY + ((0:subhologramResolution(1)-1)*slmPixelSize(1))  );                                   
                                
    % compute distance of each pixel
    r = sqrt(XKERNEL.^2 + YKERNEL.^2 + distance.^2); % Kirchhoff

    if bPhaseCompensation
        h = exp( 1i .* k .* (r-distance) ) ./ r; % this is the convolution kernel
    else
        h = exp( 1i .* k .* r ) ./ r; 
    end

    % we want this to be round
    maxValue = min(r(1,:));
    mask = double(r<=maxValue);            

    % apodize kernel
    blurSize = 20;                

    %mask = imfilter(mask, fspecial('gaussian', 3.*[blurSize blurSize], blurSize), 'same');                
    mask = conv2(mask, fspecial('gaussian', 3.*[blurSize blurSize], blurSize));                
    mask = imresize(mask, size(h));

    h = h.*mask;
    
    % for backpropagation, use conjugate transpose of kernel 
    if bBackpropagation
        h = h';
    end
end