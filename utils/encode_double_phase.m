function [dpac_3pi, dpac_2pi] = encode_double_phase(field, shift_mean_phase_DPAC,...
    max_phase, center_phase)

    if nargin < 2
		shift_mean_phase_DPAC = false;
    end
    
    if nargin < 3
        max_phase = 3 * pi;
    end

    max_field = max(abs(field(:)));
    if ( max_field > 1 )
        % if this happens, the exponents will be complex-valued, resulting in
        % phase+amplitude SLM pattern
        disp(['WARNING: image needs to be normalized for this to work! Max ' ...
        	  'value of x is ' num2str(max_field)]);
        field = field ./ max_field;
    end

    phase_target = angle(field);

    if shift_mean_phase_DPAC
        phase_target = phase_target - mean(phase_target(:)) + center_phase;
        phase_target(phase_target < -pi) = phase_target(phase_target < -pi) + 2*pi;
        phase_target(phase_target >= pi) = phase_target(phase_target >= pi) - 2*pi;
    end
    
    pa = phase_target - acos(abs(field));
    pb = phase_target + acos(abs(field));
    phi = pa;
    phi(1:2:end,2:2:end,:) = pb(1:2:end,2:2:end,:);
    phi(2:2:end,1:2:end,:) = pb(2:2:end,1:2:end,:);

    % 3 pi range (or whatever max_phase is set to)
    
    % nominally center the phases
    dpac_3pi = phi + max_phase/2;

    % wrap anything outside the 0 to max_phase range by 2pi
    while any(dpac_3pi(:) < 0)
        dpac_3pi(dpac_3pi < 0) = dpac_3pi(dpac_3pi < 0) + 2 * pi;
    end
    while any(dpac_3pi(:) >= max_phase)
        dpac_3pi(dpac_3pi >= max_phase) = ...
            dpac_3pi(dpac_3pi >= max_phase) - 2 * pi;
    end

    % normalize to [0 1]
    dpac_3pi = dpac_3pi / max_phase;   

    % 2 pi range
    if nargout > 1
        % wrap the phase values to range -pi to pi
        phi_2pi = angle(exp(1i .* phi));

        % normalize phase to range [0 1] encoding phase [0 2*pi]
        dpac_2pi = (phi_2pi + pi) ./ (2*pi);
    end

end