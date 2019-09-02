% Takes a complex field in a mat file and back propagates it 10cm to the
% SLM, followed by dpac.

clear all;
addpath('utils');

color = {'_r', '_g', '_b'};
scenes = {'landscape_0000', 'landscape_0001',...
    'landscape_day_0000', 'landscape_day_0001',...
    'landscape_night_0000', 'landscape_night_0001',...
    'lfGenLowPolyCastles_0000', 'lfGenLowPolyForest_0000'};

for scene = 1:length(scenes)

    slm_dist = 10e-2;  % 10 cm
    pixel_pitch = 6.4e-6;  % 6.4 um
    slm_resolution = [1080, 1920];
    wavelengths = [635 520 450] * 1e-9;  % rgb, 635, 520, 450 nm
    max_phases = [3 3.1 3.1] * pi;
    center_phase = pi/4;  % for better use of SLM quantization
    eyepiece_focal_length = 2.5e-2;  % 2.5 cm
    wavenums = 2*pi ./ wavelengths;

    % Loading from a light field folder
    lf_folder = ['data/' scenes{scene}];
    [lf, d] = load_unity_light_field(lf_folder, eyepiece_focal_length, [], true);

    % compute the complex field for the lf
    field = holographic_stereogram_olas(lf, d, wavenums, pixel_pitch);
    % replace the above with one of these for a different algorithm
    % % field = holographic_stereogram_hs(lf, d, wavenums);
    % % field = holographic_stereogram_pas(lf, d, wavenums);
    % % field = holographic_stereogram_apas(lf, d, wavenums);

    % for fresnel, also replace the loader so that it's rgbd input
    % % [rgb, depth] = load_unity_light_field(lf_folder, eyepiece_focal_length, [], true, true);
    % % use_fast_single_plane_fresnel_approx = false;
    % % if use_fast_single_plane_fresnel_approx
    % %     field = sqrt(rgb);
    % % else
    % %     phase_compensation_for_prop = true;
    % %     field = fresnel_hologram_rgbd(rgb, depth, wavenums, phase_compensation_for_prop, pixel_pitch);
    % % end

    % save the field too
    field_file = [lf_folder '_olas'];
    save(field_file, 'field');

    % pad up to slm resolution if it's not
    pad_amount = slm_resolution - [size(field, 1) size(field, 2)];
    if any(pad_amount > 0)
        field = padarray(field, floor(pad_amount / 2), 'pre');
        field = padarray(field, ceil(pad_amount / 2), 'post');
    end

    % propagate to slm
    field_slm = zeros(size(field));
    for c = 1:length(wavelengths)
        lambda = wavelengths(c);
        field_slm(:,:,c) = propagation_slm(field(:,:,c), [pixel_pitch pixel_pitch], ...
            lambda, -slm_dist, slm_resolution, slm_resolution);
    end

    % flip for convention
    field_slm = conj(field_slm);

    % global scaling
    field_slm = field_slm / max(field_slm(:));
    
    % dpac
    dpac_3pi = zeros(size(field_slm));
    for c = 1:length(wavelengths)
        dpac_3pi(:,:,c) = encode_double_phase(field_slm(:,:,c), true,...
            max_phases(c), center_phase);
        imwrite(dpac_3pi(:,:,c), [field_file color{c} '.png'], 'BitDepth', 8);
    end
    % save to file
    imwrite(dpac_3pi, [field_file '_rgb.png'], 'BitDepth', 8);

end
