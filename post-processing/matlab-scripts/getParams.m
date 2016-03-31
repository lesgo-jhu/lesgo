function [nx,ny,nz2,nz_tot,nproc,z_i,L_x,L_y,L_z,dx,dy,dz,u_star,nu_molec]=getParams( param_file_name )
% Reads in simulation parameters from lesgo_param.out
%   Uses regular expressions

lesgoParam = fileread(param_file_name);

tokens = regexp(lesgoParam, 'nx, ny, nz, nz_tot :\s*([^]+)', 'tokens');
params = sscanf( tokens{1}{1}, '%d' ) ;
nx = params(1);
ny = params(2); 
nz2 = params(3); 
nz_tot = params(4);

tokens = regexp(lesgoParam, 'nproc :\s*([^]+)', 'tokens');
nproc = sscanf(tokens{1}{1}, '%d' ) ;

tokens = regexp(lesgoParam, 'z_i :\s*([^]+)', 'tokens');
z_i = sscanf(tokens{1}{1}, '%f' ) ;

tokens = regexp(lesgoParam, 'L_x, L_y, L_z :\s*([^]+)', 'tokens');
params = sscanf( tokens{1}{1}, '%f' ) ;
L_x = params(1); 
L_y = params(2); 
L_z = params(3);

tokens = regexp(lesgoParam, 'dx, dy, dz :\s*([^]+)', 'tokens');
params = sscanf( tokens{1}{1}, '%f' ) ;
dx = params(1); 
dy = params(2); 
dz = params(3);

tokens = regexp(lesgoParam, 'u_star :\s*([^]+)', 'tokens');
u_star = sscanf(tokens{1}{1}, '%f' ) ;

tokens = regexp(lesgoParam, 'nu_molec :\s*([^]+)', 'tokens');
nu_molec = sscanf(tokens{1}{1}, '%f' ) ;

end

