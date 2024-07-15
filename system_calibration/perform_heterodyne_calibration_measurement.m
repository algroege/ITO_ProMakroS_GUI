function phaseResult = perform_heterodyne_calibration_measurement(app)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform_heterodyne_measurement - Brief description of the function's purpose
%
% Detailed explanation of the function and its operations, including any
% relevant background information or context.
%
% Syntax:
%   output = perform_heterodyne_measurement(app)
%
% Inputs:
%   app - reference to the assiciated ProMakroS App instance (singleton)
%
% Outputs:
%   output - Description of the output (data type, dimensions, etc.)
%   output1 - Description of the first output parameter (data type, dimensions, etc.)
%   output2 - Description of the second output parameter (data type, dimensions, etc.)
%   ...
%
% Examples:
%   Example 1:
%     % Brief description of the example
%     result = FUNCTION_NAME(input1, input2);
%
%   Example 2:
%     % Brief description of the example
%     [output1, output2] = FUNCTION_NAME(input1, input2);
%
% See also:
%   List of related functions or scripts
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

%% Methodical Parameters
N = 4; % 4-Phase Algorithm
% Fundamental wavelengths in Pixel
lam_01 = 20; % <--------- FROM APP
lam_02 = 21.581; % <--------- FROM APP
lam_03 = 23.162; % <--------- FROM APP

% lam_01 = 20; % <--------- FROM APP
% lam_02 = 40; % <--------- FROM APP
% lam_03 = 60; % <--------- FROM APP

% Schwebung 1.x
Lam_11 = lam_01 * lam_02 /(lam_02 - lam_01);
Lam_12 = lam_02 * lam_03 /(lam_03 - lam_02);
% Lam_13 = lam_01 * lam_03 /(lam_03 - lam_01);
% Schwebung 2.x
Lam_21 = Lam_11 * Lam_12 /(Lam_12 - Lam_11);
% Lam_22 = Lam_12 * Lam_13 /(Lam_13 - Lam_12);

%% Define spatial vars
x_vec = (1:app.projectorSize(1));
xy = ones(app.projectorSize(2),1) * x_vec;

%% Allocate memory for frames to capture
imagesLam_01 = zeros(app.cameraSize(2),app.cameraSize(1),N);
imagesLam_02 = imagesLam_01;
imagesLam_03 = imagesLam_01;

%% Perform Measurement
% 1. Lambda 1
phi_01_true = 2*pi*(xy)/lam_01;
waitPeriode = 0.5;  % <--------- FROM APP
% preview(app.camera)
for ii=0:N-1
    I = cos(phi_01_true-2*pi*ii/N);
    frame = uint8(rescale(I)*255);
    set(app.projectorImageHandle, 'CData', frame);
    drawnow;
    pause(waitPeriode);
    imagesLam_01(:,:,ii+1) = app.recentFrame;
end
% 2. Lambda 2
phi_02_true = 2*pi*(xy)/lam_02;
for ii=0:N-1
    I = cos(phi_02_true-2*pi*ii/N);
    frame = uint8(rescale(I)*255);
    set(app.projectorImageHandle, 'CData', frame);
    drawnow;
    pause(waitPeriode);
    imagesLam_02(:,:,ii+1) = app.recentFrame;
end
% 3. Lambda 3
phi_03_true = 2*pi*(xy)/lam_03;
for ii=0:N-1
    I = cos(phi_03_true-2*pi*ii/N);
    frame = uint8(rescale(I)*255);
    set(app.projectorImageHandle, 'CData', frame);
    drawnow;
    pause(waitPeriode);
    imagesLam_03(:,:,ii+1) = app.recentFrame;
end
% app.stopCamera(app)
% tic
%% Calculate wrapped Phase
phi_01 = atan2(imagesLam_01(:,:,4) - imagesLam_01(:,:,2), imagesLam_01(:,:,3) - imagesLam_01(:,:,1));
phi_02 = atan2(imagesLam_02(:,:,4) - imagesLam_02(:,:,2), imagesLam_02(:,:,3) - imagesLam_02(:,:,1));
phi_03 = atan2(imagesLam_03(:,:,4) - imagesLam_03(:,:,2), imagesLam_03(:,:,3) - imagesLam_03(:,:,1));
% figure(10)
% imagesc(phi_01)
% figure(20)
% imagesc(phi_02)
% figure(30)
% imagesc(phi_03)
% remove offset ?
phi_01 = phi_01 + pi;
phi_02 = phi_02 + pi;
phi_03 = phi_03 + pi;
%% Generate Modulation mask
I_n_sin_sum = imagesLam_01(:,:,1)*sin(2*pi*1/4) + imagesLam_01(:,:,2)*sin(2*pi*2/4) + imagesLam_01(:,:,3)*sin(2*pi*3/4) + imagesLam_01(:,:,4)*sin(2*pi*4/4);
I_n_cos_sum = imagesLam_01(:,:,1)*cos(2*pi*1/4) + imagesLam_01(:,:,2)*cos(2*pi*2/4) + imagesLam_01(:,:,3)*cos(2*pi*3/4) + imagesLam_01(:,:,4)*cos(2*pi*4/4);
M = 2/N*sqrt(I_n_sin_sum.^2 + I_n_cos_sum.^2);
% B_0 = imagesLam_01(:,:,2) + 1*imagesLam_01(:,:,4);
M_threshold = 10;
M_Mask = ones(size(M));
M_Mask(M < M_threshold) = NaN;

%%
% Schwebung 1.x: reines differenzsignal:
phi_11 = phi_01 - phi_02;
phi_12 = phi_02 - phi_03;
% phi_13 = phi_01 - phi_03;
% Schwebung 1.x: phasenrichtiges Schwebungssignal:
phi_11(phi_11<0) = phi_11(phi_11<0) + 2*pi;
phi_12(phi_12<0) = phi_12(phi_12<0) + 2*pi;
% phi_13(phi_13<0) = phi_13(phi_13<0) + 2*pi;
% Schwebung 1.x: Skaliere phi_1x auf gleiche Steigung wie phi_0x
% phi_11_scaled = Lam_11/lam_01 * phi_11;
% % Schwebung 1.x: Streifenordnung O_phi(x)
% O_phi_11 = round((phi_11_scaled - phi_01)/(2*pi));

% Schwebung 2.1: reines differenzsignal
phi_21 = phi_11 - phi_12;
% Schwebung 2.1: phasenrichtiges Schwebungssignal:
phi_21(phi_21<0) = phi_21(phi_21<0) + 2*pi;
% Schwebung 2.1: Skaliere phi_21 auf gleiche Steigung wie phi_11
phi_21_scaled = Lam_21/Lam_11 * phi_21;
phi_22_scaled = Lam_21/Lam_12 * phi_21;
% Schwebung 2.1: Streifenordnung O_phi(x)
O_phi_21 = round((phi_21_scaled - phi_11)/(2*pi));
O_phi_22 = round((phi_22_scaled - phi_12)/(2*pi));

% Entfalte Schwebung 1.x:
PHI_11 = phi_11 + O_phi_21 * 2*pi;
PHI_12 = phi_12 + O_phi_22 * 2*pi;
% Mittelwert
PHI_1avg = 0.5*(PHI_11 + PHI_12*Lam_12/Lam_11);

% Entfalte Phasen 0.x:
% Skaliere steigung der Schwebung 2 auf Phase
PHI_1avg_01_scaled = Lam_11/lam_01 * PHI_1avg;
PHI_1avg_02_scaled = Lam_11/lam_02 * PHI_1avg;
PHI_1avg_03_scaled = Lam_11/lam_03 * PHI_1avg;
% Streifenordnung auf basis des Mittelwerts
O_phi_11 = round((PHI_1avg_01_scaled - phi_01)/(2*pi));
O_phi_12 = round((PHI_1avg_02_scaled - phi_02)/(2*pi));
O_phi_13 = round((PHI_1avg_03_scaled - phi_03)/(2*pi));
% Entfalte Schwebung 1.x:
PHI_01 = phi_01 + O_phi_11 * 2*pi;
PHI_02 = phi_02 + O_phi_12 * 2*pi;
PHI_03 = phi_03 + O_phi_13 * 2*pi;

% Mittelwert aller Phasen
PHI_0avg = 1/3*(PHI_01 + PHI_02*lam_02/lam_01 + PHI_03*lam_03/lam_01);
% Maskieren
phaseResult = PHI_0avg .* M_Mask;
end