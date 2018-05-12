function Cs = getCSOpt2(p)
%GETCSOPT2 reads the average Cs coefficient
%   Cs = getCSOpt2(p) reads the average Cs coefficient
%   lesgo parameters are provided as struct p read using p = getParams(...)

% size of record
N = p.nx*p.ny*p.nz_tot;

% Cs coefficient
fname = './output/cs_opt2.bin';
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
Cs = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
fclose(fid);

end

