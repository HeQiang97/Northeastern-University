
function [  ] = Main(  )


clc;
close all;
clear all;

%N =100;
runs =20;
InformedAgentsSize =10;
MaximumSimulationSteps = 50;
leng = MaximumSimulationSteps + 1;
lengM=0:2:50;
mu = 0.05;

XD = 10;
UI = 2;
alpha = 1;    
Randomness = 0.1;
beta = 0.4;
% Network
% BA
% averageDegree = 12;
% A = BAnet(N,averageDegree/2,averageDegree/2);
% sp = sparse(A);
% M=sum(sum(A~=0));
% net = CreateMap(sp);

%ER  
% averageDegree = 12;
% a = averageDegree/N;
% [net,spp, A] = ErdosRenyi(N,a);
% sp = sparse(A);
% M=sum(sum(A~=0));
%  net = CreateMap(sp);
  
%WS
% averageDegree = 12;
% sp = WattsStrogatzCreator(N,averageDegree,0.1);
% A=sp;
% M=sum(sum(A~=0));
% net = CreateMap(sp);


%[a1,a2]=textread('email.txt','%d%d','headerlines',2);% good
[a1,a2]=textread('USairport500-net.txt','%d%d','headerlines',2);
%[a1,a2]=textread('facebook_combined.txt','%d%d','headerlines',2);
%[a1,a2]=textread('out-advogato.txt','%d%d','headerlines',2);%bad
%[a1,a2]=textread('NetHEPT.txt','%d%d','headerlines',2);
%[a1,a2]=textread('out.petster-friendships-hamster.txt','%d%d','headerlines',2);

 A0=[a1,a2]; 
 %A0=[a1+1,a2+1]; 
 M=length(unique(A0));
 N = max(max(A0));
 sp = sparse(A0(:,1),A0(:,2),1,N,N);
 A=full(sp); 
 MM=sum(sum(A~=0));
 net = CreateMap(sp); 

betweennesses = betweenness_centrality(sp);
[D ~] = all_shortest_paths(sp,struct('algname','floyd_warshall'));
closenesses = (N-1) ./ sum(D);
clear('sp')

meanAllOpininos0AvgMat = zeros(runs,leng);
meanAllOpininos1AvgMat = zeros(runs,leng);
meanAllOpininos2AvgMat = zeros(runs,leng);
meanAllOpininos5AvgMat = zeros(runs,leng);
followers0Avg = zeros(runs,leng);
followers1Avg = zeros(runs,leng);
followers2Avg = zeros(runs,leng);
followers5Avg = zeros(runs,leng);

for run=1:runs
    fprintf(1,'Algorithm Run %d\n',run);
    net0 = containers.Map(net.keys,net.values);
    net1 = containers.Map(net.keys,net.values);
    net2 = containers.Map(net.keys,net.values);
    net3 = containers.Map(net.keys,net.values);
    net4 = containers.Map(net.keys,net.values);
    net5 = containers.Map(net.keys,net.values);
    net6 = containers.Map(net.keys,net.values);

    x = randi(3, N, 1)-2;
    x1=sum(x==-1), x2=sum(x==0), x3=sum(x==1)
    [meanAllOpinions0,followers0] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,D,1,Randomness,alpha,beta,betweennesses,closenesses);
    [meanAllOpinions1,followers1] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,D,8,Randomness,alpha,beta,betweennesses,closenesses);
    [meanAllOpinions2,followers2] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,D,4,Randomness,alpha,beta,betweennesses,closenesses);
    [meanAllOpinions5,followers5] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,D,7,Randomness,alpha,beta,betweennesses,closenesses);
    meanAllOpininos0AvgMat(run,:) = meanAllOpinions0;
    meanAllOpininos1AvgMat(run,:) = meanAllOpinions1;
    meanAllOpininos2AvgMat(run,:) = meanAllOpinions2;
    meanAllOpininos5AvgMat(run,:) = meanAllOpinions5;
  
end
if runs ~=1 
   meanAllOpininos0AvgMat = sum(meanAllOpininos0AvgMat) / runs;
   meanAllOpininos1AvgMat= sum(meanAllOpininos1AvgMat) / runs;
   meanAllOpininos2AvgMat = sum(meanAllOpininos2AvgMat) / runs;
   meanAllOpininos5AvgMat = sum(meanAllOpininos5AvgMat) / runs;
end

end
