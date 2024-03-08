function varargout = zoomOnPlatform(viewer, platform, height)
currentPos = campos(viewer);
currentOrient = camorient(viewer);
pos = platform.Position;
campos(viewer,[pos(1) pos(2) height]);
drawnow;
[varargout{1:nargout}] = snapshot(viewer);
campos(viewer,currentPos);
camorient(viewer,currentOrient);
end