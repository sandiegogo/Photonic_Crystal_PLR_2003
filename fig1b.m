

function Main
clc;
clear;
fid = fopen('fig1b.txt','w');

n1 = 1;
n2 = 33; % unit cells's number is 16, why here is 33? why not 32? 
f1 = 1;
f2 = 50;
nf = 500;
fstep = (f2-f1)/nf;

% thickness remain the same?
dlayer(1) = 16 % 1 unit cell include 2 layers
for i = 1:1:16
    dlayer(2*i) = 4 ;
    dlayer(2*i+1) = 16 ;
end

%
for iw = 1:1:nf
    f = f1+fstep*iw;
    omega = 2*f*pi;

%materials' EM properties    
    epsilon(1) = 1;
    fmu(1) = 1;
    for i = 1:1:16
        epsilon(2*i) =  -8;
        epsilon(2*i+1) = 1;
        fmu(2*i) = -2;
        fmu(2*i+1) = 1;
    end
   
    Q = [1,0;0,1];
    
    for l=n1:1:n2-1
        [transfer] = transmtx(omega,epsilon,fmu,l);
        eps = sqrt(epsilon(l+1)) * sqrt(fmu(l+1));%can get -4 correctly, but why l+1?
        dd0 = dlayer(l+1);
        [p] = prop(fmu,omega,eps,l,dd0);%propagation matrix
        Q = Q*transfer*p;  
    end
    
     s11 = -Q(2,1)/Q(2,2);
     s21 = Q(1,1)-Q(1,2)*Q(2,1)/Q(2,2);
     t = abs(s21);% why abs?
     r = abs(s11);

s=' ';
fprintf(fid,'%f',t^2)% square means Transmittance of energy.
fprintf(fid,'%c',s)
fprintf(fid,'%f\n',f)

%draw pictures
scatter(t^2,f,'k')%why can't I draw? but I can draw the data with ORIGIN. 

end
end

function [transfer] = transmtx(omega,epsilon,fmu,l)

ckl1 = -sqrt(epsilon(l+1)/fmu(l+1));
ckl2 = -sqrt(epsilon(l)/fmu(l));
[transfer] = [ 0.50 * (1.0 + ckl2/ckl1),0.50 * (1.0 - ckl2/ckl1);0.50 * (1.0 - ckl2/ckl1), 0.50 * (1.0 + ckl2/ckl1)];
end

function [p] = prop(fmu,omega,eps,l,dd0)

% what if let c be 299.8? then f's unit will be MHz
c=299.8;
ckl = eps*omega/c; % wave number
ckl = ckl*dd0;
[p] = [exp(i*ckl),0;0,exp(-i*ckl)];
end



