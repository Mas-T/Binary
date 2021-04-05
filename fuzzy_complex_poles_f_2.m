function [f_min,f_max,alpha2]=fuzzy_complex_poles_f_2(a11,a12,T)
%find minimum and maximum values of function 
% w is the renonande frequency of the pole
% T is de maximum time to study
f_2=@(delta)(exp(-a11*T*(delta+1e-8))*sin(a12*T*(delta+1e-8)))/T/(delta+1e-8);

% the min is round 2*pi/a12
nmax=floor((T*a12-pi/2)/2/pi);
f_max=max(f_2(0),f_2(1));
f_min=min(f_2(0),f_2(1));
for n=0:nmax
   [x1,fval] = fminsearch(@(x)-f_2(x),(pi/2+2*n*pi)/a12/T);
   if (x1<1)&&(x1>0)
    f_max=max(-fval,f_max);
   end
end
% the max is round pi/a12

nmax=floor((a12*T-3*pi/2)/2/pi);

for n=0:nmax
   [x2,fval] = fminsearch(f_2,(2*n*pi+3*pi/2)/a12/T);
   if (x2<1)&&(x2>0)
     f_min=min(fval,f_min);
   end
end
alpha2=@(delta) (f_2(delta)-f_min)/(f_max-f_min); 
end