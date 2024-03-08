% Padded measurement model
function z = cvmeasPadded(x,varargin)
z = cvmeas(x,varargin{:});
numPads = 3 - size(z,1);
numStates = size(z,2);
z = [z;zeros(numPads,numStates)];
end