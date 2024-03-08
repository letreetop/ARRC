function phd = initRadarPHD(varargin)
% initcvgmphd adds 1 component when detection is an input. It adds no
% components with no inputs. 
phd = initcvgmphd(varargin{:});
phd.ProcessNoise = 5*eye(3);

% Update measurement model to use padded measurements
phd.MeasurementFcn = @cvmeasPadded;
phd.MeasurementJacobianFcn = @cvmeasjacPadded;
end