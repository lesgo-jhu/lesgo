function p = getParams(param_file_name)
% GETPARMS reads in the lesgo simulation parameters.
%   p = getParams(param_file_name) reads in simulation parameters 
%   from lesgo_param.out as a struct

lesgoParam = fileread(param_file_name);

tokens = regexp(lesgoParam, 'nx, ny, nz, nz_tot :\s*([^]+)', 'tokens');
params = sscanf( tokens{1}{1}, '%d' ) ;
p.nx = params(1);
p.ny = params(2); 
p.nz2 = params(3); 
p.nz_tot = params(4);

tokens = regexp(lesgoParam, 'nproc :\s*([^]+)', 'tokens');
p.nproc = sscanf(tokens{1}{1}, '%d' ) ;

tokens = regexp(lesgoParam, 'z_i :\s*([^]+)', 'tokens');
p.z_i = sscanf(tokens{1}{1}, '%f' ) ;

tokens = regexp(lesgoParam, 'L_x, L_y, L_z :\s*([^]+)', 'tokens');
params = sscanf( tokens{1}{1}, '%f' ) ;
p.L_x = params(1); 
p.L_y = params(2); 
p.L_z = params(3);

tokens = regexp(lesgoParam, 'dx, dy, dz :\s*([^]+)', 'tokens');
params = sscanf( tokens{1}{1}, '%f' ) ;
p.dx = params(1); 
p.dy = params(2); 
p.dz = params(3);

tokens = regexp(lesgoParam, 'u_star :\s*([^]+)', 'tokens');
p.u_star = sscanf(tokens{1}{1}, '%f' ) ;

tokens = regexp(lesgoParam, 'nu_molec :\s*([^]+)', 'tokens');
p.nu_molec = sscanf(tokens{1}{1}, '%f' ) ;

tokens = regexp(lesgoParam, 'write_endian :\s*([^]+)\n', 'tokens');
tc = textscan(tokens{1}{1},'%s',1);
p.write_endian = char(tc{1});


tokens = regexp(lesgoParam, 'read_endian :\s*([^]+)\n', 'tokens');
tc = textscan(tokens{1}{1},'%s',1);
p.read_endian = char(tc{1});

% Set the format specifier for fread
if strcmp(p.write_endian,'BIG_ENDIAN')
    p.fmt = 's';
else
    p.fmt = 'a';
end

% Set the indices for each processor
dummynproc = 1:p.nproc;
p.zmin_buf = dummynproc*(p.nz2-1)-p.nz2+2;
p.zmax_buf = dummynproc*(p.nz2-1)+1;

% build grid
p.x    = 0.0    : p.dx : p.L_x-p.dx+1E-6;
p.y    = 0.0    : p.dy : p.L_y-p.dy+1E-6;
p.z_w  = 0.0    : p.dz : p.L_z+1E-6;           % for avg vels, avg (rs's, dudz,dvdz,txz,tyz)
p.z_uv = p.dz/2 : p.dz : p.L_z+p.dz/2+1E-6;     % for inst vels, avg (txx,txy,tzz,txy)

end

