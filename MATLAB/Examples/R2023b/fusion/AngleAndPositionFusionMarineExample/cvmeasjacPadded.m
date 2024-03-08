% Jacobian of padded measurement model
function H = cvmeasjacPadded(x,varargin)
H = cvmeasjac(x,varargin{:});
numPads = 3 - size(H,1);
numStateElements = size(H,2);
H = [H;zeros(numPads,numStateElements)];
end