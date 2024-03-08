function updateDisplay(viewer, scenario, detections, tracks)

% Define position and velocity selectors
posSelector = [1 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 1 0];
velSelector = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1];

% Color order
colorOrder = [1.0000    1.0000    0.0667
    0.0745    0.6235    1.0000
    1.0000    0.4118    0.1608
    0.3922    0.8314    0.0745
    0.7176    0.2745    1.0000
    0.0588    1.0000    1.0000
    1.0000    0.0745    0.6510];

% Plot Platforms
plotPlatform(viewer, [scenario.Platforms{:}],'ECEF',Color=colorOrder(2,:),TrajectoryMode='Full');

% Plot Coverage
covConfig = coverageConfig(scenario);
covConfig(2).Range = 10e3;
plotCoverage(viewer,covConfig,'ECEF',Color=colorOrder([1 7],:),Alpha=0.2);

% Plot Detections
sIdx = cellfun(@(x)x.SensorIndex,detections);
plotDetection(viewer, detections(sIdx == 1), 'ECEF',Color=colorOrder(1,:));
helperPlotAngleOnlyDetection(viewer, detections(sIdx == 2),Color=colorOrder(7,:));

% Plot Tracks
plotTrack(viewer,tracks,'ECEF',PositionSelector=posSelector,...
    VelocitySelector=velSelector,LineWidth=2,Color=colorOrder(4,:));
drawnow;

end