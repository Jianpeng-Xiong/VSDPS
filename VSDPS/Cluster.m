%% The clustering method used on dominated solutions
function [SPoint,Cflag,C] = Cluster(Global,Population,I,Pre_population)
         K1=Global.M+1;
         K2=Global.M*3;
         K=K1+floor(I*(K2-K1));
         if K>K2
            K=K2;
         end
         Pop_Decs=Population.decs;
         Pre_Pop_Decs=Pre_population.decs;
         C=zeros(K,Global.D);
         indexs=Farthest_first_selection(Population,K);
         SPoint=Pop_Decs(indexs,:);
         Dis=pdist2(SPoint,Pop_Decs);
         [~,Cflag]=min(Dis,[],1);
         for i=1:K
             SPoint(i,:)=mean(Pop_Decs(Cflag==i,:));
         end
         for i=1:K
             C(i,:)=SPoint(i,:)-Pre_Pop_Decs(indexs(i),:);
         end
end

