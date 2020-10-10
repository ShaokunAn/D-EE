clear all
s = RandStream('mcg16807','Seed',29); RandStream.setGlobalStream(s);
Y=csvread('hspcs_pca.csv',0,0);
N=size(Y,1);
d=2;
X0 = 1e-5*randn(N,d);
csvwrite("hspcs_matlab_X0.csv",X0);