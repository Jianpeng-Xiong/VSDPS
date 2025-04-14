function VSDPS(Global)
% <algorithm> <V>
% A dynamic multi-objective evolutionary algorithm with variable stepsize and dual prediction strategies
% nr    ---   2 --- Maximum number of solutions replaced by each offspring
% LArchive --- 150 --- Maximum of the Cur_Archive

% This is a simple demo of VSDPS
%--------------------------------------------------------------------------------------------------------
% If you find this code useful in your work, please cite the following paper "Hu Penga, Chen Pi, Jianpeng 
% Xiong, Debin Fan, Fanfan Shen. A dynamic multi-objective evolutionary algorithm with variable stepsize and 
% dual prediction strategies. Future Generation Computer Systems, 2024".
%--------------------------------------------------------------------------------------------------------
% This function is implmented by Chen Pi
%--------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
% More information can visit Hu Peng's homepage: https://whuph.github.io/index.html
%--------------------------------------------------------------------------------------------------------
%------------------------------- DMO --------------------------------------------------------------------

  %% Parameter setting
   [delta,nr,LArchive,proM,disM] = Global.ParameterSet(0.9,2,150,4,2);

  %% Generate the weight vectors
   [W,Global.N] = UniformPoint(Global.N,Global.M);
   T = 30;     
   
  %% Detect the neighbours of each solution
   B = pdist2(W,W);
   [~,B] = sort(B,2);
   B = B(:,1:T);
   
  %% Generate random population
   Global.Initialization();
   Z = min(Global.Population.objs,[],1);
   
  %% Optimization
    % The historical nondominated solutions is stored in the Cur_Archive
      Cur_Archive=[];
      while Global.NotTermination() 
            % Detecting changes in the environment using re-evaluation methods
             if Global.hasChanged()  
                   % Cur_non is non-dominated solutions of the current population
                     [FrontNo,~] = NDSort(Global.Population.objs,Global.Population.cons,Global.N);
                     Cur_non=Global.Population(FrontNo==1);
             if Global.gen-50==Global.problem.tauT   
                   % The first environmental change occurs, no response is made, and the population is re-evaluated.
                     Pre_population=Global.Population;
                     C = mean(Global.Population.decs);          
                     Global.Population = First_to_reaction(Global, Global.Population,W,Z,proM,disM);
             elseif Global.gen-50>Global.problem.tauT 
                     Cur_population=Global.Population;
                   % Variable stepsize
                     I=Intense_detect(Global,Pre_population,Global.Population);
                     C = mean(Global.Population.decs);       
                     [SPoint,Cflag,Center]=Cluster(Global,Global.Population,I,Pre_population); 
                   % Dual prediction strategies
                     Global.Population = DualPredictionStrategies(Global,Global.Population,pre_C,C,W,Z,Cur_Archive,Cflag,Center,SPoint);
                     Pre_population=Cur_population;
             end
                     Cur_Archive= [Cur_Archive,Cur_non];
             if length(Cur_Archive) > LArchive
                     Cur_Archive= Cur_Archive(:,length(Cur_Archive)-LArchive+1:end);
             end
                 pre_C  = C;  
                 Z = min(Global.Population.objs,[],1);   
             end
        % improved MOEA/D-DE
          for i = 1 : Global.N
              if rand < delta
                 PN = B(i,randperm(end));
              else
                 PN = randperm(Global.N);
              end
             % Generate an offspring
              if rand>0.4
                     Offspring = DE(Global.Population(i),Global.Population(PN(1)),Global.Population(PN(2)));
              else
              if mod(i,3)==0
                     Offspring=DE_current_to_lbest(Global.Population,Global.Population(i),W,Z,PN,Global.lower,Global.upper);
              elseif mod(i,3)==1
                     Offspring=GAhalf([Global.Population(i) Global.Population(PN(1))]);
              elseif mod(i,3)==2
                     Offspring = DE_2best(Global.Population,Global.Population(i),W,Z,PN,Global.lower,Global.upper);
              end
              end
             % Update the ideal point
               Z = min(Z,Offspring.obj);
             % Update the solutions in P by Tchebycheff approach
               g_old = max(abs(Global.Population(PN).objs-repmat(Z,length(PN),1)).*W(PN,:),[],2);
               g_new = max(repmat(abs(Offspring.obj-Z),length(PN),1).*W(PN,:),[],2);
               Global.Population(PN(find(g_old>=g_new,nr))) = Offspring;
          end 
      end 
   end
          