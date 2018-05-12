function [u,v,w] = getAvgVelUV(p)
%GETAVGVELUV reads the average velocity in the uv grid
%   [u,v,w] = getAvgVelUV(p) reads the average velocity on the uv grid
%   lesgo parameters are provided as struct p read using p = getParams(...)

% size of record
N = p.nx*p.ny*p.nz_tot;

% Velocity
fname = './output/veluv_avg.bin';
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
u = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
v = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
w = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
fclose(fid);

end

