function  [meanMajOpinions,followers,x1] = MySimulationBodyMethodMulti(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,alpha,targets,R)
N=size(x,1);
x(targets)=1;
meanMajOpinions=[];
followers=[];
H=zeros(N,N);
T = zeros(N,1);
for it=1:MaximumSimulationSteps
    H(:,it)=x;
    cnt = 0;
    cnp = 0;
    cn0 = 0;
    TT=3;
    Q = T;
    for i=1:N
        if (i== targets)
           continue; 
        end
        lastX = x;             
        adj = find(A(i,:)==1);
        if length(adj) == 0
           continue;
        end
        index = adj(randi(length(adj)));
        Na = zeros(N,1);
        Ar = zeros(N,1);
        O=[+1, 0, -1]; 
        reward =zeros(N,1);     
        for j=1:3      
            if(O(j) - x(i) == 0)
              reward(j)=1;
           end
           if(O(j) - x(i) ~= 0)
              reward(j)=0;
           end
           if(O(j) == lastX(index))
             reward(j) =  reward(j)+1;
           end
           if(O(j) - lastX(index) ~= 0)
             reward(j) =  reward(j)-1;
           end 
           Q(j) = Q(j) + mu * (reward(j) - Q(j));    % error
        end
        [max_num,max_index]=max(exp(Q/TT)./sum(exp(Q(j)/TT)));
        TT=TT-0.05;
        if (TT<1)
           TT=1;
        end  
        T(i)=Q(max_index);
        if (max_index==1)
            x(i)=+1;
            cnt =cnt+1;
        end
        if (max_index==2)
           x(i)=0;
           cn0 =cn0+1;
        end
        if (max_index==3)
            x(i)=-1;
            cnp =cnp+1;
        end
        Q = zeros(N,1); 
       if (it <= MaximumSimulationSteps)
          for j=1:3
             for t=1:it
                if (O(j)==H(i,t))
                   Na(j) = Na(j)+1;  
                end
             end
             if (Na(j) == 0)
                Ar(j)=0;
             end
             if (Na(j) ~= 0)
                Ar(j) = Na(j);
             end
          end
       end
      [max_nu,max_outdex]=max(Ar);  
      if (max_outdex == max_index)
          mu=0.01;
      end
      if (max_outdex ~= max_index)
          mu=0.05;
      end                            
    end
    cnf=cnt-cnp;
    meanMajOpininos(it) = cnf;
    followersratio(it) = cnt/N;
    meanMajOpinions=[meanMajOpinions; meanMajOpininos(it)];
    followers=[followers; followersratio(it)];
    x1=x;
end
end
