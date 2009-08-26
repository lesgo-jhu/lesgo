clear all;
close all;

ngen=1;

if ngen == 1 
  nfiles = 5;
  tstart = 0;
  fname = ['check_ke.out.200k';
         'check_ke.out.400k';
         'check_ke.out.600k';
         'check_ke.out.800k';
         'check_ke.out.1000k']
elseif ngen == 3  
nfiles = 10;
tstart =  1000000*2e-4;
fname = ['check_ke.out.50k';
         'check_ke.out.60k';
         'check_ke.out.70k';
         'check_ke.out.80k';
         'check_ke.out.90k';
         'check_ke.out.100k';
         'check_ke.out.110k';
         'check_ke.out.120k';
         'check_ke.out.130k';
         'check_ke.out.140k']
end

for nf=1:nfiles
  raw_data(:,:,nf) = load(strtrim(fname(nf,:)));
end

isize = size(raw_data,1);
time = zeros(isize*nfiles,1);
ke = zeros(isize*nfiles,1);


for nf=1:nfiles
  ioffset = (nf - 1)*isize;
  for i=1:isize
    time(i+ioffset) = raw_data(i,1,nf);
    ke(i+ioffset) = raw_data(i,2,nf);
  end
end

time = time + tstart;

plot(time,ke)
%hold on;
%xlabel('time')
%ylabel('KE')
%title('checkbd-3 : ngen1 : 600k iterations');
%print -dpng ke_checkbd-3_ngen1_600k.png
