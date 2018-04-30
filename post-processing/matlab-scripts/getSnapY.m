function [u,v,w] = getSnapY(p,step,loc)
%GETSNAPY reads y-plane snapshot
%   getSnapY(p,step) y-plane snapshot at timestep step and y location loc.
%   Lesgo parameters are provided as struct p read using p = getParams(...)

% size of record
N = p.ny*p.nz_tot;

% Velocity
fname = ['./output/vel.y-',sprintf('%0.5f',loc),'.',num2str(step),'.bin'];
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
u = reshape(fread(fid,N,'double',p.fmt),p.ny,p.nz_tot);
v = reshape(fread(fid,N,'double',p.fmt),p.ny,p.nz_tot);
w = reshape(fread(fid,N,'double',p.fmt),p.ny,p.nz_tot);
fclose(fid);

end

