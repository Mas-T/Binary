function [f_min,f_max,alpha2]=fuzzy_complex_poles_f_1(a11,a12,T)
%find minimum and maximum values of function 
% w is the renonande frequency of the pole
% T is de maximum time to study
f_1=@(delta)(exp(-a11*T*(delta+1e-8))*cos(a12*T*(delta+1e-8))-1)/T/(delta+1e-8);

% the min is round 2*pi/a12
nmax=floor(T*a12/2/pi);
f_max=max(f_1(0),f_1(1));
f_min=min(f_1(0),f_1(1));
for n=0:nmax
   [x1,fval] = fminsearch(@(x)-f_1(x),2*n*pi/a12/T);
   if (x1<1)&&(x1>0)
    f_max=max(-fval,f_max);
   end
end
% the max is round pi/a12

nmax=floor((a12*T-pi)/2/pi);

for n=0:nmax
   [x2,fval] = fminsearch(f_1,(2*n*pi+pi)/a12/T);
   if (x2<1)&&(x2>0)
     f_min=min(fval,f_min);
   end
end
alpha2=@(delta) (f_1(delta)-f_min)/(f_max-f_min);
end