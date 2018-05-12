function [uu,vv,ww,uw,vw,uv] = getReyStress(p)
%GETREYSTRESS reads the average stress
%   [uu,vv,ww,uw,vw,uv] = getReyStress(p) reads the Reynolds stress from the 
%   time averaged data. lesgo parameters are provided as struct p read using 
%    p = getParams(...)

% size of record
N = p.nx*p.ny*p.nz_tot;

% Reynolds stress
fname = './output/rs.bin';
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
uu = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
vv = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
ww = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
uw = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
vw = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
uv = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
fclose(fid);
    
end

