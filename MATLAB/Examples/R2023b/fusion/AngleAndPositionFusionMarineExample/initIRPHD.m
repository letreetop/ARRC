function phd = initIRPHD(varargin)
% This adds no components
phd = initcvgmphd; 
phd.ProcessNoise = 5*eye(3);

% Update measurement model to use padded measurements
phd.MeasurementFcn = @cvmeasPadded; 
phd.MeasurementJacobianFcn = @cvmeasjacPadded;
end