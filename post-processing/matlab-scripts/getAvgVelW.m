function w = getAvgVelW(p)
%GETAVGVELW reads the average velocity in the w grid
%   w = getAvgVelW(p) reads the average velocity on the w grid
%   lesgo parameters are provided as struct p read using p = getParams(...)

% size of record
N = p.nx*p.ny*p.nz_tot;

% Velocity
fname = './output/velw_avg.bin';
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
w = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
fclose(fid);

end

