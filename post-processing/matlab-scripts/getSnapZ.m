function [u,v,w] = getSnapZ(p,step,loc)
%GETSNAPZ reads z-plane snapshot
%   getSnapZ(p,step) z-plane snapshot at timestep step and z location loc.
%   Lesgo parameters are provided as struct p read using p = getParams(...)

% size of record
N = p.nx*p.ny;

% Velocity
fname = ['./output/vel.z-',sprintf('%0.5f',loc),'.',num2str(step),'.bin'];
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
u = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny);
v = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny);
w = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny);
fclose(fid);

end

