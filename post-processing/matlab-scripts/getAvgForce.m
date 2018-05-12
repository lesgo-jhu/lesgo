function [fx,fy,fz] = getAvgForce(p)
%GETAVGFORCE reads the average force
%   [fx,fy,fz] = getAvgForce(p) reads the average force. 
%   lesgo parameters are provided as struct p read using p = getParams(...)

% size of record
N = p.nx*p.ny*p.nz_tot;

% Force
fname = './output/force_avg.bin';
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
fx = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
fy = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
fz = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
fclose(fid);

end

