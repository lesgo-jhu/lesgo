function [u,v,w] = getSnapX(p,step,loc)
%GETSNAPX reads x-plane snapshot
%   getSnapX(p,step) x-plane snapshot at timestep step and x location loc.
%   Lesgo parameters are provided as struct p read using p = getParams(...)

% size of record
N = p.nx*p.nz_tot;

% Velocity
fname = ['./output/vel.x-',sprintf('%0.5f',loc),'.',num2str(step),'.bin'];
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
u = reshape(fread(fid,N,'double',p.fmt),p.nx,p.nz_tot);
v = reshape(fread(fid,N,'double',p.fmt),p.nx,p.nz_tot);
w = reshape(fread(fid,N,'double',p.fmt),p.nx,p.nz_tot);
fclose(fid);

end

