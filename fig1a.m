clear;
clc;
fid = fopen('fig_1_a_data.txt','w')
f=linspace(0,50,1000);%f is a vector
%K is pure imaginary, its real part is zero
K1=imag((acos(1+1/4.*(sin(8*pi.*f/75)).^2)))/pi;
%because function'y=cos(Ka)'is an even function.
K2=-K1;
plot(K1,f,'-');
hold on;
plot(K2,f,'-'); 
hold on;
freal=[0,75/8,75/4,225/8,75/2,375/8];
Kreal=[0,0,0,0,0,0];
plot(Kreal,freal,'*')
axis([-1,1,0,50])

%wrting arrays into a file.
s=' ';
[b1 b2]=size(K1);
for i=1:b1
    for j=1:b2
       fprintf(fid,'%10d',f(i,j));
       fprintf(fid,'%c',s);
       fprintf(fid,'%10d',K1(i,j));
       fprintf(fid,'%c',s);
       fprintf(fid,'%10d\n',K2(i,j));  

    end
end
