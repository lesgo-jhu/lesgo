function [u2,v2,w2,uw,vw,uv] = getAvgVel2(p)
%GETAVGTAU reads the average squared velocites
%   [u2,v2,w2,uw,vw,uv] = getAvgVel2(p) reads the average squared velocities
%   lesgo parameters are provided as struct p read using p = getParams(...)

% size of record
N = p.nx*p.ny*p.nz_tot;

% Velocity^2
fname = './output/vel2_avg.bin';
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
u2 = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
v2 = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
w2 = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
uw = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
vw = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
uv = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
fclose(fid);

end

