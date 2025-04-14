function flag = Crowding(PopObj,FrontNo)
% Calculate the crowding distance of each solution front by front
    [N,M]    = size(PopObj);
    flag=zeros(1,N);
    Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts)
        if sum(flag)==0.8*N
                break;
        end
        Front = find(FrontNo==Fronts(f));
        for i = 1 : M
            [~,Rank] = sortrows(PopObj(Front,i));
             if sum(flag)==0.8*N
                break;
            end
            flag(Front(Rank(1)))=1;
             if sum(flag)==0.8*N
                break;
            end
            flag(Front(Rank(1)))=1;
             if sum(flag)==0.8*N
                break;
            end
            for j = 2 : length(Front)-1
                flag(Front(Rank(j)))=1;
                 if sum(flag)==0.8*N
                break;
                end
            end
        end
    end
end