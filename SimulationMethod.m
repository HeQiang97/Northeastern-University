
function [meanMajOpinions,followers] = SimulationMethod( x,InformedAgentsSize,MaximumSimulationSteps,mu,A,D,mode,Randomness,alpha,beta,betweennesses,closenesses)

N = size(x,1);
XD = 1;
meanMajO=[];
follow=[];
indegs = sum(A);
outdegs = sum(A');
R=2;
Steps=1;
cnt = 0;
cnp = 0;
cn0 = 0;
meanMajOpinions =zeros(length(InformedAgentsSize)+1,1);
for i=1:(N-InformedAgentsSize)
    if (x(i)==1)
	    cnt =cnt+1;
    end 
    if (x(i)==0)
		cn0 =cn0+1;
    end 
    if (x(i)==-1)
	    cnp =cnp+1;
    end                        
end    
 
cnf=cnt-cnp;
majorityMeanOpinions0 = cnf;
initFol=cnt/N;
sp = sparse(A);

if mode == 1
   P=zeros(N,1);
   targets=[];
   for r=1:R
       for i = 1 : N
         SP=zeros(N,1);
         if outdegs(i)==0
            continue;
         end
         AOUT1 = find(A(i,:)==1);
         Leht=length(find(A(i,:)==1));
         Ain = find(A(:,i)==1);
         for j = 1 : Leht
            Aout= find(A(AOUT1(j),:)==1);
            B=intersect(Ain,Aout);
            C=union(Ain,Aout);
            D=(length(B)/length(C));
            if D==0 
               SP(i)=SP(i)+x(j)+x(i)*outdegs(i)/(outdegs(i)+outdegs( AOUT1(j)));
            end
            if D~=0
               SP(i)=SP(i)+x(j)+x(i)*D;
            end
         end
         P(i)=x(i)+SP(i); 
       end
       [OB11,OB1] = sort(P,'descend');
       targets = OB1(1:ceil(InformedAgentsSize/R))';
       x(targets)=XD;
       [meanMajOpinionss,followerss,x1] = MySimulationBodyMethodMulti1(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,alpha,targets,R);
       x=x1;
   end
   meanMajO =[meanMajO;meanMajOpinionss];
   follow= [follow;followerss];
   meanMajOpinions=[majorityMeanOpinions0;meanMajO];
   followers = [initFol;follow];
end
if mode == 2
   EPN2_measure = zeros(N,1);
   eps = 0.0001;
   for i = 1 : N       
      f1 = ((outdegs(i)) / (indegs(i) + outdegs(i)))
      ss = 0;   
      for j = 1 : outdegs(i)
          s(j) = 1/2* outdegs1(j) +   1/6*betweennesses(j) +  1/3*closenesses(j);
          ss= ss +  ( beta * s(j) + (1-beta) * mean(s(A(j,:) == 1)) )      
      end
      f2 = ss ./ outdegs(i)
      f2 = f2 ./ max(s);
      EPN2_measure(i) = f1 * f2;
   end
   [~,EP2] = sort(EPN2_measure,'descend');
   targets = EP2(1:InformedAgentsSize)'; %  
   x(targets)=XD;
   [meanMajOpinions,followers] = MySimulationBodyMethodh(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,alpha,targets);
   meanMajOpinions = [majorityMeanOpinions0;meanMajOpinions];
   followers = [initFol;followers];
end
   
if mode == 3 
     EPN1_measure = zeros(N,1);
     eps = 0.0001;
     for i = 1 : N              
         s(i) = 1/2* outdegs1(i) +   1/6*betweennesses(i) +  1/3*closenesses(i);
         f1 = ((outdegs(i) + eps) / (indegs(i) + outdegs(i)));
         f2 = (mean(s(A(i,:) == 1)) / max(indegs));
         EPN1_measure(i) = f1 * f2;
     end
     [~,EP1] = sort(EPN1_measure,'descend');
     targets = EP1(1:InformedAgentsSize)';
     x(targets)=XD;
     [meanMajOpinions,followers] = MySimulationBodyMethodh(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,alpha,targets);
     meanMajOpinions = [majorityMeanOpinions0;meanMajOpinions];
    followers = [initFol;followers];
end

if mode == 4
   [~,idxx] = sort(indegs,'descend');
   targets = idxx(1:InformedAgentsSize)';
   x(targets)=XD;
   [meanMajOpinionss,followerss] = MySimulationBodyMethodh(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,alpha,targets);       
   meanMajO =[meanMajO;meanMajOpinionss];
   follow= [follow;followerss];
   meanMajOpinions=[majorityMeanOpinions0;meanMajO];
   followers = [initFol;follow];  
end
      
if mode == 5
   randdata=randperm(N);
   targets=randdata(1:InformedAgentsSize(agent));
   x(targets)=XD;            
   [meanMajOpinions,followers] = MySimulationBodyMethodh(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,alpha,targets);
   meanMajOpinions = [majorityMeanOpinions0;meanMajOpinions];
   followers = [initFol;followers];
end

if mode == 6
   [~,xx] = sort(closenesses,'descend');
   targets = xx(1:InformedAgentsSize)';
   x(targets)=XD;
end

if mode == 7
   M=zeros(1,N);
   M1=zeros(N,N);
   M2=zeros(N,N);
   CC=zeros(N,1);
   P=zeros(N,N);
   for i=1:N
       M(i)=sum(A(i,:))+sum(A(:,i));
       for j=1:N
           P(i,j)=(A(i,j)+A(j,i))./ M(i);
           for q=1:N
               M1(i,j)=M1(i,j)+ P(i,q).*P(q,j); 
           end
           M2(i)=M2(i)+(P(i,j)+M1(i,j))^2;  
       end
       if indegs(i) ~=0
          CC(i)=  M2(i)./indegs(i);  
       end       
    end 
    [~,USED] = sort(CC,'descend');
    targets = USED(1:InformedAgentsSize)';
    x(targets)=XD;  
    [meanMajOpinionss,followerss] = MySimulationBodyMethodh(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,alpha,targets);       
    meanMajO =[meanMajO;meanMajOpinionss];
    follow= [follow;followerss];
    meanMajOpinions=[majorityMeanOpinions0;meanMajO];
    followers = [initFol;follow];    
end
if mode == 8
   for i= 1: N
       t=0.5;
       for j = 1:N
           Page(j) = (A(i,j).*D(i,j))/outdegs(j);
        end
        PageRank(i)= (1-t)/N + t*Page(i);
   end
   [~,xxxx] = sort(PageRank,'descend');
   targets = xxxx(1:InformedAgentsSize)';
   x(targets)=XD;                     
   [meanMajOpinionss,followerss] = MySimulationBodyMethodh(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,alpha,targets);
   meanMajO =[meanMajO;meanMajOpinionss];
   follow= [follow;followerss];
   meanMajOpinions=[majorityMeanOpinions0;meanMajO];
   followers = [initFol;follow];     
end
end
               
