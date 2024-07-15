function pts_w = backward_core(K, E, pts_i)
    % Author: Alex G.
    % Description: core backward projection mode for a single image with
    % distortion.    
    % Inputs:
    %   K - Intrinsic matrix
    %   E - Extrinsic matrix
    %   v_distort - distortion coefficient vector
    %   pts_w - world coordinates of points
    %   distort_model - 'none', 'rad_tan' or 'full' (default)
    %   distort_plane - 'normal' (defualt) or 'pixel'
    % Outputs:
    %   pts_i - image coordinates of points. N x 2
    
    % Forward projection
    N = size(pts_i,1);

    % 1. transform from pixel to camera coordinate system
    pts_i_temp = [pts_i ones(N,1)];
    pts_pre = pts_i_temp/K;
    % 2. apply distortion
    pts_distort = pts_pre;
    % 3. Transform to world coordinate
    pts_temp = pts_distort/E;
    pts_w = pts_temp;           
end