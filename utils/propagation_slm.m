%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Implementation of free space propagation of coherent light with a 
%   phase-only spatial light modulator.
%
%   Input parameters (all distances in meters!):
%
%       fieldIn         - complex-valued field on the input plane
%       slmPixelSize    - [height width] of an SLM pixel
%       slmResolution   - [height width] in number of pixels of the SLM
%       imageResolution - width and height of target image (should be square)
%       lambda          - wavelength of light
%       distance        - propagation distance (use positive value for
%                           propagation from SLM to target plane and negative
%                           value for backpropagation)
%       
%   Output parameters:
%       fieldOut        - complex-valued field on the output plane
%
%   Other parameters:
%       bIsFresnel      - use Fresnel propagation instead of Kirchhoff-Fresnel
%
%   Gordon Wetzstein
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fieldOut = propagation_slm(fieldIn, slmPixelSize, lambda, distance, slmResolution, imageResolution)
        
    % check if we're propagating back to SLM plane
    bBackpropagation = true;
    if distance>=0
        bBackpropagation = false;
    end
    distance = abs(distance);
        
    % wave number
    k = 2*pi/lambda;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % compute size of subhologram via grating equation 
    
    % max diffraction angle 
    maxDiffractionAngleRad = asin(lambda ./ (2*slmPixelSize));
    
    % size of subhologram in m
    subhologramSize = 2 * distance * tan(maxDiffractionAngleRad);
    
    %subhologramSize = 2 * distance * lambda ./ (2*pixelSize); % this is the small angle approximation
    
    % subhologram feature size in m  
    subhologramFeatureSize = subhologramSize ./ imageResolution;
    
    % resolution of subhologram in number of slm pixels
    subhologramResolution = round((subhologramSize ./ slmPixelSize));
    
    % this should be even
    subhologramResolution = floor(subhologramResolution./2) .* 2;
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % compute convolution kernel
    
    % get coordinates of pixel centers in kernel
    startY = -slmPixelSize(1)*(subhologramResolution(1)-1)/2;
    startX = -slmPixelSize(2)*(subhologramResolution(2)-1)/2;
    [XKERNEL, YKERNEL] = meshgrid(   startX + ((0:subhologramResolution(2)-1)*slmPixelSize(2)), ...
                                    startY + ((0:subhologramResolution(1)-1)*slmPixelSize(1))  );                                   
               
    % compute distance of each pixel 
    r = sqrt(XKERNEL.^2 + YKERNEL.^2 + distance.^2); % Kirchhoff        


    h = exp( 1i .* k .* (r-distance) ) ./ r; % this is the convolution kernel

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
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fieldOut = propagation_conv(fieldIn, h);
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function fieldOut = propagation_conv(fieldIn, h)       

    opFFT2  = @(x) fft2(x);
    opIFFT2 = @(x) ifft2(x);

    % pad input image by a bit more than half the kernel size
    padSize = floor(size(h)./2)+1;
    fieldIn = padarray(fieldIn,padSize, median(abs(fieldIn(:))));                        

    % get kernel in frequency domain
    H = psf2otf(h, size(fieldIn));       

    % get OTF and perform convolution in frequency domain
    fieldOut = opIFFT2( opFFT2(fieldIn).*H );      

    % crop back
    if prod(padSize)>0
        fieldOut = fieldOut( padSize(1)+1:end-padSize(1), padSize(2)+1:end-padSize(2) );
    end   

end