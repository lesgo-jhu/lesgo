clear;
close all;

nproc=8;
zbottom=0.31746;
sets=[200:200:4000]';
nsets=length(sets);
for ns=1:nsets
  istart=0;
  iend=0;
  for n=1:nproc
    fname=strcat("uvw_avg_z-",int2str(sets(1)),"k-",int2str(sets(ns)),"k.dat.c",int2str(n-1));

    fid=fopen(fname,'r');

    dat = fscanf(fid, '%g %g %g %g ', [4 inf]);
    fclose(fid);
    dat=dat';
    ndat=size(dat,1);
    istart=iend+1;
    iend = istart + ndat - 1 -1; % additional -1 to skip redundant value a top of domain slice
    z(ns,istart:iend) = dat(1:ndat-1,1);
    u(ns,istart:iend) = dat(1:ndat-1,2);
    v(ns,istart:iend) = dat(1:ndat-1,3);
    w(ns,istart:iend) = dat(1:ndat-1,4);
  end
  
end
npoints=size(z,2);
for ns=1:nsets
  for n=1:npoints
    if(u(nsets,n) == 0.)
      u_err(ns,n)=0.;
    else
      u_err(ns,n) = (u(ns,n) - u(nsets,n))/u(nsets,n);
    end

    if(v(nsets,n) == 0.)
      v_err(ns,n)=0.;
    else
      v_err(ns,n) = (v(ns,n) - v(nsets,n))/v(nsets,n);
    end

    if(w(nsets,n) == 0.)
      w_err(ns,n)=0.;
    else
      w_err(ns,n) = (w(ns,n) - w(nsets,n))/w(nsets,n);
    end
  end
end

%  Remove no flow bottom region
z=z-zbottom;

for ns=1:nsets
  fname=strcat("uvw_avg_z-",int2str(sets(1)),"k-",int2str(sets(ns)),"k.dat");
  fid=fopen(fname,'w');
  for n=1:npoints
    if(z(ns,n) >= 0)
      fprintf(fid,'%g \t %g \t %g \t %g \t %g \t %g \t %g\n',z(ns,n), u(ns,n), v(ns,n), w(ns,n), u_err(ns,n), v_err(ns,n), w_err(ns,n));
    end
  end
  fclose(fid);  
end
