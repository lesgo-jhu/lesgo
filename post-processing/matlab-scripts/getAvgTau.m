function [txx,txy,tyy,txz,tyz,tzz] = getAvgTau(p)
%GETAVGTAU reads the average stress
%   [txx,txy,tyy,txz,tyz,tzz] = getAvgTau(p) reads the average stress
%   lesgo parameters are provided as struct p read using p = getParams(...)

% size of record
N = p.nx*p.ny*p.nz_tot;

% Stress
fname = './output/tau_avg.bin';
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
txx = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
txy = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
tyy = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
txz = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
tyz = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
tzz = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
fclose(fid);

end

