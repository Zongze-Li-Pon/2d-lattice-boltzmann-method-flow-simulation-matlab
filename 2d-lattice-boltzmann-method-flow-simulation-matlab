%% ========================================================================
%  2D Lattice Boltzmann Method Flow Simulation in MATLAB (D2Q9, BGK)
%  Karman Vortex Street Simulation
%  
%  This script demonstrates a complete 2D lattice Boltzmann simulation
%  for incompressible channel flow with:
%    - pressure inlet / pressure outlet (Zou/He)
%    - no-slip top and bottom walls with bounce-back boundary
%    - circular obstruction with bounce-back boundary
%
%  The purpose of this script is educational:
%  with only a few hundred lines of MATLAB code, one can run a 2D fluid
%  simulation and visualize velocity/pressure fields in real time.
%
%  Author: Zongze Li
%  ========================================================================

%% Clean workspace
clear, close all, clc;

%% ========================================================================
%  1. Domain and flow parameters
%  ========================================================================

% Lattice/domain size
lx = 400;                      % number of lattice nodes in x direction
ly = 200;                      % number of lattice nodes in y direction

% Flow condition
Re   = 200;                    % Reynolds number
uMax = 0.15;                   % target maximum inlet velocity

% Circular obstruction geometry
obst_x = lx / 4;               % x location of cylinder center
obst_y = ly / 2;               % y location of cylinder center
obst_d = ly / 10;              % cylinder diameter


%% ========================================================================
%  2. LBM parameters
%  ========================================================================
 
% Kinematic viscosity in lattice unit
% Here the characteristic length is the cylinder diameter
nu = uMax * obst_d / Re;

% Relaxation time and relaxation parameter
tau   = 3 * nu + 0.5;
omega = 1 / tau;

fprintf('Relaxation time tau = %.6f\n', tau);

% Stability recommended working range check
if tau > 1.5 || tau < 0.51
    fprintf('\nTau is out of the recommended range: 0.51 < tau < 1.5\n');
    return;
end

% Simulation length and visualization frequency
maxT  = 800000;                % total time steps
tPlot = 5;                    % plot every tPlot steps


%% ========================================================================
%  3. D2Q9 lattice constants
%  ========================================================================

% Lattice weights
w = [4/9, ...
     1/9, 1/9, 1/9, 1/9, ...
     1/36, 1/36, 1/36, 1/36];

% Discrete lattice velocities
cx = [ 0,  1,  0, -1,  0,  1, -1, -1,  1];
cy = [ 0,  0,  1,  0, -1,  1,  1, -1, -1];

% Opposite directions (used in bounce-back)
opp = [1, 4, 5, 2, 3, 8, 9, 6, 7];


%% ========================================================================
%  4. Initial distribution function
%  ========================================================================

% Initialize distribution function with equilibrium at rho = 1, u = 0
fIn = reshape(w' * ones(1, lx * ly), 9, lx, ly);

% Allocate equilibrium and post-collision distributions
fEq  = zeros(9, lx, ly);
fOut = zeros(9, lx, ly);


%% ========================================================================
%  5. Analytical Poiseuille profile (for reference/comparison)
%  ========================================================================

% Interior y locations used for analytical inlet profile. Because
% bounce-back wall presents at half way of the lattice link, so the
% "actual" wall location is at y = 0.5 and ly-0.5.
yLine = 0.5 : ly-0.5;
H = ly - 1;

% Parabolic inlet profile (used only as a theoretical comparison curve)
uInPro = 4 * uMax / (H^2) * ((yLine - 0.5) .* H - (yLine-0.5).^2);


%% ========================================================================
%  6. Pressure difference and density boundary values
%  ========================================================================

% For plane Poiseuille flow:
dP = 8 * uMax * lx * nu / ly^2;

% In LBM, p = cs^2 * rho = (1/3) * rho  -->  dP = (1/3) * dRho
dRho = 3 * dP;

% Set inlet/outlet density
rho_in  = 1.0 + 0.5 * dRho;
rho_out = 1.0 - 0.5 * dRho;


%% ========================================================================
%  7. Convergence monitor storage
%  ========================================================================

% Preallocate convergence history
iter      = zeros(maxT, 1);
ConvergeX = zeros(maxT, 1);
ConvergeY = zeros(maxT, 1);


%% ========================================================================
%  8. Geometry and plotting coordinates
%  ========================================================================

% Mesh coordinates for plotting and geometry definition
[yGrid, xGrid] = meshgrid(1:ly, 1:lx);

% Circular obstruction mask
obst = (xGrid - obst_x).^2 + (yGrid - obst_y).^2 <= (obst_d / 2)^2;

% Linear indices of obstruction nodes (used for bounce-back and masking)
bbRegion = find(obst);

% Linear indices of bottom and top walls (used for bounce-back)
bottomWall = false(lx, ly);
bottomWall(:,1) = true;
bbBottom = find(bottomWall);

topWall = false(lx, ly);
topWall(:,ly) = true;
bbTop = find(topWall);

%% ========================================================================
%  9. Initial macroscopic field
%  ========================================================================

% Start from a Poiseuille-like velocity field across the domain
% This helps the simulation begin closer to the expected flow state.
Hinit   = ly - 1;
y_phys  = yGrid - 0.5;

ux = 4 * uMax / (Hinit^2) * (y_phys .* Hinit - y_phys .* y_phys);
uy = zeros(lx, ly);
rho = 1.0;

% Initialize zero velocity on all bounce-back boundaries.
ux(bbRegion) = 0;
uy(bbRegion) = 0;
ux(bbBottom) = 0;
uy(bbBottom) = 0;
ux(bbTop) = 0;
uy(bbTop) = 0;

% Build initial equilibrium distribution from (rho, ux, uy)
for i = 1:9
    cu = 3 * (cx(i) * ux + cy(i) * uy);
    fIn(i,:,:) = rho .* w(i) .* ...
        (1 + cu + 0.5 * (cu.^2) - 1.5 * (ux.^2 + uy.^2));
end


%% ========================================================================
%  10. Initialize visualization figures
%  ========================================================================

% Set global figure location and size.
set(groot,'defaultFigurePosition',[100 100 600 500]);

% ---- Figure 1: convergence history ----
fig1 = figure("WindowStyle","normal");
ax1 = axes(fig1);
hConvX = semilogy(ax1, nan, nan, '.k');
hold(ax1, 'on');
hConvY = semilogy(ax1, nan, nan, '.g');
hold(ax1, 'off');
title(ax1, 'Convergence Check');
xlabel(ax1, 'Iteration');
ylabel(ax1, 'Relative L2 Change');
legend(ax1, 'u_x', 'u_y');

% ---- Figure 2: velocity magnitude contour ----
fig2 = figure("WindowStyle","normal");
ax2 = axes(fig2);
uMag0 = zeros(lx, ly);
uMag0(bbRegion) = nan;
hImgU = imagesc(ax2, uMag0');
title(ax2, 'Velocity Contour');
axis(ax2, 'equal');
axis(ax2, 'off');
colorbar(ax2,'southoutside');
axis(ax2,'tight')
clim(ax2,[0 0.6*uMax]);

% ---- Figure 3: velocity vectors ----
fig3 = figure("WindowStyle","normal");
ax3 = axes(fig3);
skipY = round(ly/15);
skipX = skipY;
hQuiver = quiver(ax3, ...
    xGrid(1:skipX:lx,1:skipY:ly), ...
    yGrid(1:skipX:lx,1:skipY:ly), ...
    zeros(size(xGrid(1:skipX:lx,1:skipY:ly))), ...
    zeros(size(yGrid(1:skipX:lx,1:skipY:ly))));
hold(ax3, 'on');
plot(ax3, [1,lx], [1,1], 'k', [1,lx], [ly,ly], 'k', 'LineWidth', 2);

% Plot obstruction outline
theta = linspace(0, 2*pi, 200);
xObst = obst_x + (obst_d/2) * cos(theta);
yObst = obst_y + (obst_d/2) * sin(theta);
plot(ax3, xObst, yObst, 'r', 'LineWidth', 2);

hold(ax3, 'off');
title(ax3, 'Velocity Vector');
xlim(ax3, [1, lx]);
ylim(ax3, [1, ly]);
axis(ax3, 'equal');

% ---- Figure 4: compare profiles with analytical Poiseuille curve ----
fig4 = figure("WindowStyle","normal");
ax41 = subplot(2,2,1, 'Parent', fig4);
hP11 = plot(ax41, nan, nan, 'k', 'LineWidth', 2);
hold(ax41, 'on');
hP12 = plot(ax41, nan, nan, 'ob', 'MarkerSize', 3);
hold(ax41, 'off');
title(ax41, 'x = 1');
axis(ax41, 'equal');

ax42 = subplot(2,2,2, 'Parent', fig4);
hP21 = plot(ax42, nan, nan, 'k', 'LineWidth', 2);
hold(ax42, 'on');
hP22 = plot(ax42, nan, nan, 'ob', 'MarkerSize', 3);
hold(ax42, 'off');
title(ax42, sprintf('x = obstruction (x/L = %.3f)', obst_x/lx));
axis(ax42, 'equal');

ax43 = subplot(2,2,3, 'Parent', fig4);
hP31 = plot(ax43, nan, nan, 'k', 'LineWidth', 2);
hold(ax43, 'on');
hP32 = plot(ax43, nan, nan, 'ob', 'MarkerSize', 3);
hold(ax43, 'off');
title(ax43, 'x = lx / 2');
axis(ax43, 'equal');

ax44 = subplot(2,2,4, 'Parent', fig4);
hP41 = plot(ax44, nan, nan, 'k', 'LineWidth', 2);
hold(ax44, 'on');
hP42 = plot(ax44, nan, nan, 'ob', 'MarkerSize', 3);
hold(ax44, 'off');
title(ax44, 'x = lx');
axis(ax44, 'equal');

% ---- Figure 5: pressure contour ----
fig5 = figure("WindowStyle","normal");
ax5 = axes(fig5);
pField0 = zeros(lx, ly);
pField0(bbRegion) = nan;
hImgP = imagesc(ax5, pField0');
axis(ax5, 'equal');
title(ax5, 'Pressure Contour');
colorbar(ax5);
pMean  = 1/3;
pRange = 1.5*dRho;
clim(ax5,[pMean-pRange pMean+pRange]);

% ---- Figure 6: centerline pressure ----
fig6 = figure("WindowStyle","normal");
ax6 = axes(fig6);
hCenterP = plot(ax6, nan, nan, 'LineWidth', 1.5);
hold(ax6, 'on');

% Mark obstruction projection on the centerline
xObstLeft  = (obst_x - obst_d/2) / lx;
xObstRight = (obst_x + obst_d/2) / lx;
xline(ax6, xObstLeft,  '--r', 'LineWidth', 1.5);
xline(ax6, xObstRight, '--r', 'LineWidth', 1.5);

hold(ax6, 'off');
title(ax6, 'Centerline Pressure');
xlabel(ax6, 'x / L');
ylabel(ax6, 'p / p_{in}');
grid(ax6, 'on');

% ---- Figure 7: vorticity contour ----
fig7 = figure("WindowStyle","normal");
ax7 = axes(fig7);

vort0 = zeros(lx,ly);
vort0(bbRegion) = nan;

hImgVort = imagesc(ax7, vort0');
axis(ax7,'equal')
axis(ax7,'off')
title(ax7,'Vorticity Contour')
colorbar(ax7,'southoutside');
axis(ax7,'tight')
clim(ax7,[-uMax/obst_d uMax/obst_d])
colormap(ax7,turbo)


%% ========================================================================
%  10.5 GIF output settings
%  ========================================================================

saveGif = true;                 % true to save GIFs
gifStartCycle = 1;              % can change to e.g. 2000 if you only want later flow
gifStep = 10;                   % save one GIF frame every gifStep visualization updates
gifDelay = 0.05;                % GIF playback delay time (seconds)

velGifFile  = 'velocity_contour.gif';
vortGifFile = 'vorticity_contour.gif';

velGifInitialized  = false;
vortGifInitialized = false;


%% ========================================================================
%  11. Main time-stepping loop
%  ========================================================================

for cycle = 1:maxT

    %% --------------------------------------------------------------------
    %  11.1 Recover macroscopic variables from distribution functions
    %  --------------------------------------------------------------------
    rho = sum(fIn, 1);
    ux  = reshape(cx * reshape(fIn, 9, lx * ly), 1, lx, ly) ./ rho;
    uy  = reshape(cy * reshape(fIn, 9, lx * ly), 1, lx, ly) ./ rho;

    % Save "before boundary update" velocity for convergence check
    uxA = reshape(ux, lx, ly);
    uyA = reshape(uy, lx, ly);


    %% --------------------------------------------------------------------
    %  11.2 Collision step (BGK)
    %  --------------------------------------------------------------------
    for i = 1:9
        cu = 3 * (cx(i) * ux + cy(i) * uy);

        fEq(i,:,:) = rho .* w(i) .* ...
            (1 + cu + 0.5 * (cu.^2) - 1.5 * (ux.^2 + uy.^2));

        fOut(i,:,:) = fIn(i,:,:) - omega .* (fIn(i,:,:) - fEq(i,:,:));
    end


    %% --------------------------------------------------------------------
    %  11.3 Bounce-back on circular obstruction
    %  --------------------------------------------------------------------
    for i = 1:9
        % Bounce-back on circular obstruction
        fOut(i, bbRegion) = fIn(opp(i), bbRegion);

        % Bounce-back on bottom wall
        fOut(i, bbBottom) = fIn(opp(i), bbBottom);

        % Bounce-back on top wall
        fOut(i, bbTop) = fIn(opp(i), bbTop);
    end


    %% --------------------------------------------------------------------
    %  11.4 Streaming step
    %  --------------------------------------------------------------------
    for i = 1:9
        fIn(i,:,:) = circshift(fOut(i,:,:), [0, cx(i), cy(i)]);
    end


    %% --------------------------------------------------------------------
    %  11.5 Pressure inlet boundary (Zou/He)
    %  Left boundary: prescribe density, reconstruct unknown populations
    %  according to Zou/He's paper.
    %  --------------------------------------------------------------------
    rho(:,1,2:ly-1) = rho_in;

    ux(:,1,2:ly-1) = 1 - ...
        (sum(fIn([1,3,5],1,2:ly-1),1) + 2 * sum(fIn([4,7,8],1,2:ly-1),1)) ...
        ./ rho(:,1,2:ly-1);

    uy(:,1,2:ly-1) = 0;

    fIn(2,1,2:ly-1) = fIn(4,1,2:ly-1) + 2/3 * rho(:,1,2:ly-1) .* ux(:,1,2:ly-1);

    fIn(6,1,2:ly-1) = fIn(8,1,2:ly-1) ...
        + 1/2 * (fIn(5,1,2:ly-1) - fIn(3,1,2:ly-1)) ...
        + 1/2 * rho(:,1,2:ly-1) .* uy(:,1,2:ly-1) ...
        + 1/6 * rho(:,1,2:ly-1) .* ux(:,1,2:ly-1);

    fIn(9,1,2:ly-1) = fIn(7,1,2:ly-1) ...
        + 1/2 * (fIn(3,1,2:ly-1) - fIn(5,1,2:ly-1)) ...
        - 1/2 * rho(:,1,2:ly-1) .* uy(:,1,2:ly-1) ...
        + 1/6 * rho(:,1,2:ly-1) .* ux(:,1,2:ly-1);


    %% --------------------------------------------------------------------
    %  11.6 Pressure outlet boundary (Zou/He)
    %  Right boundary: prescribe density, reconstruct unknown populations 
    %  according to Zou/He's paper.
    %  --------------------------------------------------------------------
    rho(:,lx,2:ly-1) = rho_out;

    ux(:,lx,2:ly-1) = -1 + ...
        (sum(fIn([1,3,5],lx,2:ly-1),1) + 2 * sum(fIn([2,6,9],lx,2:ly-1),1)) ...
        ./ rho(:,lx,2:ly-1);

    uy(:,lx,2:ly-1) = 0;

    fIn(4,lx,2:ly-1) = fIn(2,lx,2:ly-1) - 2/3 * rho(:,lx,2:ly-1) .* ux(:,lx,2:ly-1);

    fIn(8,lx,2:ly-1) = fIn(6,lx,2:ly-1) ...
        + 1/2 * (fIn(3,lx,2:ly-1) - fIn(5,lx,2:ly-1)) ...
        - 1/2 * rho(:,lx,2:ly-1) .* uy(:,lx,2:ly-1) ...
        - 1/6 * rho(:,lx,2:ly-1) .* ux(:,lx,2:ly-1);

    fIn(7,lx,2:ly-1) = fIn(9,lx,2:ly-1) ...
        + 1/2 * (fIn(5,lx,2:ly-1) - fIn(3,lx,2:ly-1)) ...
        + 1/2 * rho(:,lx,2:ly-1) .* uy(:,lx,2:ly-1) ...
        - 1/6 * rho(:,lx,2:ly-1) .* ux(:,lx,2:ly-1);


    %% --------------------------------------------------------------------
    %  11.7 Recover macroscopic variables again after BC update
    %  --------------------------------------------------------------------
    rho = sum(fIn, 1);
    ux  = reshape(cx * reshape(fIn, 9, lx * ly), 1, lx, ly) ./ rho;
    uy  = reshape(cy * reshape(fIn, 9, lx * ly), 1, lx, ly) ./ rho;

    uxB = reshape(ux, lx, ly);
    uyB = reshape(uy, lx, ly);


    %% --------------------------------------------------------------------
    %  11.8 Convergence check using L2 relative change
    %  --------------------------------------------------------------------
    numeratorx   = (uxB - uxA).^2;
    denominatorx = uxA.^2;

    numeratory   = (uyB - uyA).^2;
    denominatory = uyA.^2;

    ConvergexCheck = sqrt(sum(numeratorx(:)) / sum(denominatorx(:)));
    ConvergeyCheck = sqrt(sum(numeratory(:)) / sum(denominatory(:)));

    iter(cycle)      = cycle;
    ConvergeX(cycle) = ConvergexCheck;
    ConvergeY(cycle) = ConvergeyCheck;


    %% --------------------------------------------------------------------
    %  11.9 Visualization
    %  --------------------------------------------------------------------
    if mod(cycle, tPlot) == 1

        % ---- Figure 1: convergence history ----
        set(hConvX, 'XData', iter(1:cycle), 'YData', ConvergeX(1:cycle));
        set(hConvY, 'XData', iter(1:cycle), 'YData', ConvergeY(1:cycle));

        % ---- Figure 2: velocity magnitude contour ----
        uMag = sqrt(uxB.^2 + uyB.^2);
        uMag(bbRegion) = nan;
        set(hImgU, 'CData', uMag');

        % ---- Figure 3: velocity vectors ----
        set(hQuiver, ...
            'UData', uxB(1:skipX:lx,1:skipY:ly), ...
            'VData', uyB(1:skipX:lx,1:skipY:ly));

        % ---- Figure 4: compare profiles with analytical Poiseuille curve ----
        set(hP11, 'XData', uxB(1,:) / uMax,             'YData', (1:ly)/ly);
        set(hP12, 'XData', uInPro / uMax,               'YData', (1:ly)/ly);

        set(hP21, 'XData', uxB(round(obst_x),:) / uMax, 'YData', (1:ly)/ly);
        set(hP22, 'XData', uInPro / uMax,               'YData', (1:ly)/ly);

        set(hP31, 'XData', uxB(round(lx/2),:) / uMax,   'YData', (1:ly)/ly);
        set(hP32, 'XData', uInPro / uMax,               'YData', (1:ly)/ly);

        set(hP41, 'XData', uxB(lx,:) / uMax,            'YData', (1:ly)/ly);
        set(hP42, 'XData', uInPro / uMax,               'YData', (1:ly)/ly);

        % ---- Figure 5: pressure contour ----
        rhoField  = reshape(sum(fIn,1), lx, ly);
        pField    = (1/3) * rhoField;
        pField(bbRegion) = nan;
        set(hImgP, 'CData', pField');

        % ---- Figure 6: centerline pressure ----
        set(hCenterP, ...
            'XData', (1:lx)/lx, ...
            'YData', pField(:, ly/2) / pField(1, ly/2));
        
        % ---- Figure 7: vorticity contour ----
        dudy = diff(uxB,1,2);
        dvdx = diff(uyB,1,1);

        omegaV = zeros(lx,ly);
        omegaV(1:end-1,1:end-1) = dvdx(:,1:end-1) - dudy(1:end-1,:);

        omegaV(bbRegion) = nan;

        set(hImgVort,'CData',omegaV');

        drawnow limitrate;

        
        % ---- Save GIFs: velocity contour and vorticity contour ----
        if saveGif && cycle >= gifStartCycle && mod(cycle, gifStep*tPlot) == 1

            % Save velocity contour GIF
            frameVel = getframe(fig2);
            imgVel = frame2im(frameVel);
            [imindVel, cmVel] = rgb2ind(imgVel, 256);

            if ~velGifInitialized
                imwrite(imindVel, cmVel, velGifFile, 'gif', ...
                    'Loopcount', inf, 'DelayTime', gifDelay);
                velGifInitialized = true;
            else
                imwrite(imindVel, cmVel, velGifFile, 'gif', ...
                    'WriteMode', 'append', 'DelayTime', gifDelay);
            end

            % Save vorticity contour GIF
            frameVort = getframe(fig7);
            imgVort = frame2im(frameVort);
            [imindVort, cmVort] = rgb2ind(imgVort, 256);

            if ~vortGifInitialized
                imwrite(imindVort, cmVort, vortGifFile, 'gif', ...
                    'Loopcount', inf, 'DelayTime', gifDelay);
                vortGifInitialized = true;
            else
                imwrite(imindVort, cmVort, vortGifFile, 'gif', ...
                    'WriteMode', 'append', 'DelayTime', gifDelay);
            end

        end
    end

end
