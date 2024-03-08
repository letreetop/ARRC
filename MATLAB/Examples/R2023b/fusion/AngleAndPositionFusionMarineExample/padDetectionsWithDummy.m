function detections = padDetectionsWithDummy(detections)
for i = 1:numel(detections)
    % Infrared detections have SensorIndex = 2
    if detections{i}.SensorIndex == 2
        % Add dummy measurement
        detections{i}.Measurement = [detections{i}.Measurement;0];

        % A covariance value of 1/(2*pi) on the dummy measurement produces
        % the same Gaussian likelihood as the original measurement without
        % dummy
        detections{i}.MeasurementNoise = blkdiag(detections{i}.MeasurementNoise,1/(2*pi));
    end
end