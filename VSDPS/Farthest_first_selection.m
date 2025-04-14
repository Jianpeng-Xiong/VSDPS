%% Furthest-first-choice method: Selection of better-distributed non-dominated solutions
function indexs=Farthest_first_selection(Population,t_size)
   remain=false(1,length(Population));
   [~,min_index]=min(Population.objs,[],1);
   remain(min_index)=true;
   temp_pop=Population(remain);
   while sum(remain)<t_size
       p=find(remain==false);
       Dist=pdist2(Population(p).objs,temp_pop.objs);
       [Sort_Dist,~]=sort(Dist,2);
       temp_Dist=Sort_Dist(:,1);
       [~,ins]=max(temp_Dist);
       select_index=p(ins);
       remain(select_index)=true;
       temp_pop=[temp_pop,Population(select_index)];   
   end
   indexs=find(remain==true);
end