% file "plot_powspe.m"
N=65536;
dt=0.001; %[sec]
df=1/(N*dt) % [Hz]
fmax=1/(2*dt) %[Hz]
f=[df:df:(fmax-df)]; % [Hz]
w=f*2*pi();
filename='data/2013-02-01_16-09__0.dat';
temp=importdata(filename,'',6);
S=temp.data(1:(length(temp.data())-1));
plot(f(1:500),S(1:500),'b-')
