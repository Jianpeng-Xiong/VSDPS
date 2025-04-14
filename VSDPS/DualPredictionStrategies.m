%% Dual Prediction Strategies
   function Population = DualPredictionStrategies(Global,Population,pre_C,post_C,W1,Z,Archive,Cflag,Center,SPoint)
       W = 0.5;
       OffVel=zeros(Global.N,Global.D);
       r1  = repmat(rand(Global.N,1),1,Global.D);
       r2  = repmat(rand(Global.N,1),1,Global.D);
       Pop_Decs=Population.decs;
       [FrontNo,~] = NDSort(Population.objs,Population.cons,Global.N);
     % Tchebycheff values are calculated for each individual in the population and sorted in ascending order
       TchValue=max(abs(Population.objs-repmat(Z,length(Population),1)).*W1,[],2);
       [Sort_TchValue,~]=sort(TchValue,1,'ascend');
       halfN = ceil(Global.N*0.6);
       V=post_C - pre_C;
       [Front,~] = NDSort(Archive.objs,Archive.cons,length(Archive));
       DA_index=find(Front==1);
    % The selection of the Gbest: Non-dominated individuals in the archive are selected as Gbest according to crowding distance from largest to smallest
       [Front,~] = NDSort(Archive(Front==1).objs,Archive(Front==1).cons,length(Archive((Front==1))));
       Crowdis=CrowdingDistance(Archive(Front==1).objs,Front);
       [~,Cur_index]=sort(Crowdis,2,'descend');
       P=Cur_index(:,1:floor(end/2));
       for i = 1:Global.N
           if TchValue(i) < Sort_TchValue(halfN)    
            % The top 60% of individuals use improved linear prediction strategy
              if FrontNo(i)==1
                 Pop_Decs(i,:)= Pop_Decs(i,:) +V ;
              else
                 Pop_Decs(i,:)= Pop_Decs(i,:) +Center(Cflag(i),:);
              end
           else      
         % The left 40% of individuals use dynamic particle swarm prediction strategy
           if length(P)*0.7<1
              k=1;
           else
              k=randi(floor(length(P)*0.7));
           end
         % The selection of the Pbest
           Distance=pdist2(Pop_Decs(i,:),SPoint);
           [~,j]=min(Distance,[],2);
           [~,DF]=find(Cflag==j&FrontNo==1);
           if ~isempty(DF) 
              [Front_DF,~] = NDSort(Population(DF).objs,Population(DF).cons,length(DF));
              Crowdis=CrowdingDistance(Population(DF).objs,Front_DF);
              [~,index]=sort(Crowdis,2,'descend');
              if length(DF)*0.5<1
                 m=1;
              else
                 m=randi(floor(length(DF)*0.5));
              end
              if FrontNo(i)==1
                 OffVel(i,:)= W.*V+r1(i,:).*(Pop_Decs(DF(index(m)),:)+V-Pop_Decs(i,:))+r2(i,:).*(Archive(DA_index(Cur_index(P(k)))).decs+V-Pop_Decs(i,:));
              else
                 OffVel(i,:)= W.*V+r1(i,:).*(Pop_Decs(DF(index(m)),:)+Center(Cflag(i),:)-Pop_Decs(i,:))+r2(i,:).*(Archive(DA_index(Cur_index(P(k)))).decs+V-Pop_Decs(i,:));
              end
          else 
               [~,DK]=find(Cflag==j);
               if ~isempty(DK)
                  if length(DK)==1
                     m=1;
                  else
                     m=randi(length(DK));
                  end
                  if FrontNo(i)==1
                     OffVel(i,:)= W.*V+r1(i,:).*(Pop_Decs(DK(m),:)+V-Pop_Decs(i,:))+r2(i,:).*(Archive(DA_index(Cur_index(P(k)))).decs+V-Pop_Decs(i,:));
                  else
                     OffVel(i,:)= W.*V+r1(i,:).*(Pop_Decs(DK(m),:)+Center(Cflag(i),:)-Pop_Decs(i,:))+r2(i,:).*(Archive(DA_index(Cur_index(P(k)))).decs+V-Pop_Decs(i,:));
                  end
              else 
                  index=find(FrontNo==1);
                  L_A=length(index);
                  m=randi(L_A);
                  if FrontNo(i)==1
                     OffVel(i,:)= W.*V+r1(i,:).*(Pop_Decs(index(m),:)+V-Pop_Decs(i,:))+r2(i,:).*(Archive(DA_index(Cur_index(P(k)))).decs+V-Pop_Decs(i,:));
                  else
                     OffVel(i,:)= W.*V+r1(i,:).*(Pop_Decs(index(m),:)+Center(Cflag(i),:)-Pop_Decs(i,:))+r2(i,:).*(Archive(DA_index(Cur_index(P(k)))).decs+V-Pop_Decs(i,:));
                  end
               end
           end
           Pop_Decs(i,:) = Pop_Decs(i,:) + 0.5*OffVel(i,:);
           end
       end
       Pop_Decs=Boundary_Repair(Pop_Decs,Global.lower,Global.upper);
       Population=INDIVIDUAL(Pop_Decs,OffVel);
    end