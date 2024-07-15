function telecError = evaluate_telecentricity_error(app, folderName, dlg)
    %% 1. Load Checkerboar Images
    images = read_images(folderName,'.bmp');
    
    %% 2. find points
    % Find key Points
    % patch_width = 3.0; % Kantenlänge in mm
    [imagePoints, boardSize] = detectCheckerboardPoints(images{1});
    num_pts = size(imagePoints,1);
    
    % Get board geometry
    img_num = length(images);
    batch_img_pts = zeros(num_pts, 2, img_num);
    batch_img_pts(:,:,1) = imagePoints;
    
    for k=2:img_num
        dlg.Value = k/img_num;
        % Find checkerboard points
        [imagePoints,~] = detectCheckerboardPoints(images{k});
        batch_img_pts(:,:,k) = imagePoints; % grid points in one image
    end
    
    %% 3. plot points as lines in 3d
    imagesc(images{1}, Parent=app.PhaseCalibrationAxesC1);
    hold(app.PhaseCalibrationAxesC1, "on");
    for k=1:img_num
        scatter(app.PhaseCalibrationAxesC1, batch_img_pts(:,1,k), batch_img_pts(:,2,k), 50)
    end
    hold(app.PhaseCalibrationAxesC1, "off");
    %% 4 find opposite corners
    % compute diagonal distance of corners
    distance = NaN(img_num,1);
    for k=1:img_num
        [point1 idx1] = min(batch_img_pts(:,:,k));
        [point2 idx2] = max(batch_img_pts(:,:,k));
        distance(k) = sqrt((point2(1) - point1(1))^2 + (point2(2) - point1(2))^2);
    end
    % Calculate the reference distance
    reference_distance = mean(distance);
    
    % Compute the absolute error for each element
    absolute_errors = abs(distance - reference_distance);
    
    telecError = max(absolute_errors)/reference_distance*100;
end
