function linFit = phase_to_height_calibration(calibrationMeasurementData, z_vec, dlg)
    % perform_heterodyne_measurement - Brief description of the function's purpose
    %
    % Detailed explanation of the function and its operations, including any
    % relevant background information or context.
    %
    % Syntax:
    %   output = perform_heterodyne_measurement(app)
    %
    % Inputs:
    %   calibrationMeasurementData - ...
    %   z_vec - ...
    %
    % Outputs:
    %   linFit - Matrix(m x n x k) with linear fit coefficients for each pixel
    %
    %
    % Authors:
    %   Alexander Gröger (contact@amg-optics.de)
    %
    % Version History:
    %   v1.0 - Initial version (30.06.2024)
    %
    % License:
    %   Brief description of the licensing terms or reference to a LICENSE file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calc linear polynom mx+b for each pixel: faster
    linFit = NaN(size(calibrationMeasurementData, 1), size(calibrationMeasurementData, 2), 2);
    n = length(z_vec);
    for k=1:size(calibrationMeasurementData, 1)
        dlg.Value = k/size(calibrationMeasurementData, 1);
        for j=1:size(calibrationMeasurementData, 2)
            pixel_vec = squeeze(calibrationMeasurementData(k,j,:));
            % Lineare Regression durchführen: Designmatrix erstellen
            X = [pixel_vec(:), ones(n, 1)];
            % Lineare Regression mit Backslash-Operator
            coeffs = X \ z_vec(:);
            % Koeffizienten der gefitteten Gerade
            linFit(k,j,1) = coeffs(1); % Steigung
            linFit(k,j,2) = coeffs(2); % y-Achsenabschnitt
        end
    end
    % optional: Filter data
    sizeA1 = [size(calibrationMeasurementData, 1) size(calibrationMeasurementData, 2)];
    linFit_c1 = linFit(:,:,1);
    linFit_c1_vec = linFit_c1(:);
    linFit_c1_vec = filloutliers(linFit_c1_vec, "nearest");
    linFit_c1_orem = reshape(linFit_c1_vec, sizeA1);
    linFit(:,:,1) = linFit_c1_orem;
end