% file "plot_powspe.m"

clear all

filename='data/2013-02-12_17-38__1.dat';
temp=importdata(filename,'\t',1);

dt=temp.data(1,1);
N=temp.data(1,2);
df=1.0/(N*dt); % [Hz]
fmax=1.0/(2.0*dt); %[Hz]
f=[df:df:fmax]; % [Hz]
w=f*2.0*pi(); % [rad/s]

temp=importdata(filename,'',7);
S=temp.data();

% plotting original data
figure
plot(f(1:10000),S(1:10000),'b-')

% computing averaged data for given binsize
binsize=10;
f2=[(binsize/2*df):(binsize*df):(fmax-binsize/2*df)];
S2=zeros(1,length(f2));
std_dev_S2=zeros(1,length(f2));

for i=0:(length(f2)-1)
	tmp=S((i*binsize+1):(i+1)*binsize);
	S2(i+1)=mean(tmp);
	std_dev_S2(i+1)=std(tmp);
end

% plotting averaged data
figure
plot(f2(1:1000),S2(1:1000),'k-')

