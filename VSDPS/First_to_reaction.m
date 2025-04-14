%% Respond to the first changes in the environment
function Population = First_to_reaction(Global,Population,W,Z,proM,disM)
    TchValue=max(abs(Population.objs-repmat(Z,length(Population),1)).*W,[],2);
    [Sort_TchValue,~]=sort(TchValue);
    halfN = ceil(Global.N/2);
    for i = 1:Global.N
        if TchValue(i) < Sort_TchValue(halfN) 
            Population(i) = INDIVIDUAL(Population(i).dec);
        else
            Population(i)=PNM(Population(i),{proM,disM});
%             newOneDec=Global.lower+rand(1,Global.D).*(Global.upper-Global.lower);
%             Population(i)=INDIVIDUAL(newOneDec);
        end
    end
end

