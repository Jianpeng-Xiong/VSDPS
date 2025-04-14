%% DE_more
function Offspring = DE_more(Parent,Parameter)
    %% Parameter setting
    if nargin > 3
        [CR,F,proM,disM] = deal(Parameter{:});
    else
        [CR,F,proM,disM] = deal(1,0.5,1,20);
    end
    if isa(Parent(1),'INDIVIDUAL') 
        calObj  = true;            
        parent_size=length(Parent);
        p=nchoosek(1:parent_size,3);
        parents_number=randperm(size(p,1),parent_size);
        p=p(parents_number,:);
        Parent1 = Parent(p(:,1)).decs;
        Parent2 = Parent(p(:,2)).decs;
        Parent3 = Parent(p(:,3)).decs;
    else
        calObj = false;       
    end
    [N,D]  = size(Parent1);      
    Global = GLOBAL.GetObj();

    %% Differental evolution 
    Site = rand(N,D) < CR;    
    Offspring       = Parent1;
    Offspring(Site) = Offspring(Site) + F*(Parent2(Site)-Parent3(Site)); 

    %% Polynomial mutation  
    Lower = repmat(Global.lower,N,1);
    Upper = repmat(Global.upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    if calObj
        Offspring = INDIVIDUAL(Offspring);
    end
end