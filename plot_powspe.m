% file "plot_powspe.m"
filename='data/2013-02-08_15-05__0.dat';
temp=importdata(filename,'\t',1);

dt=temp.data(1,1);
N=temp.data(1,2);
df=1.0/(N*dt); % [Hz]
fmax=1.0/(2.0*dt); %[Hz]
f=[df:df:(fmax-df)]; % [Hz]
w=f*2.0*pi();

temp=importdata(filename,'',7);
S=temp.data(1:(length(temp.data())-1));
plot(f(1:10000),S(1:10000),'b-')
