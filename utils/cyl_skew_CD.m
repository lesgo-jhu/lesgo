%clear all;
%close all;

d = 28.8*4./185.;
l = 50.4*4./185.;
ntrunk=3;
ngen=1;
scale_fact=0.5;
Ap = scale_fact^(2*(ngen-1))*d*l*ntrunk^ngen;
nfiles = 2;

fname(1,:) = ['cyl_skew_CD.out.g1.c0'];
fname(2,:) = ['cyl_skew_CD.out.g1.c1'];

for n=1:nfiles
  data(:,:,n) = load(fname(n,:));
end

nsize=size(data,1);
fD=zeros(nsize,1);
CD=zeros(nsize,1);
%  Add forces; CD
for n=1:nfiles
  fD = fD + data(:,3,n);
  CD = CD + data(:,2,n);
end
CD = CD/Ap;

fD_mean=sum(fD)/nsize
CD_mean=sum(CD)/nsize

fid=fopen('cyl_skew_CD_mean.dat.g1.c0','w');
fprintf(fid,'%s\t%i\t%f\t%f\t%f\n','Nsamples,Ap,fD,CD :',nsize,Ap,fD_mean,CD_mean);
fclose(fid);


