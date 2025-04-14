%% Intense_detect
function I = Intense_detect(Global,Pre_population,Population)
         Pop_objs=Population.objs;
         Pre_Pop_objs=Pre_population.objs;
         Min=min(Pop_objs);
         Max=max(Pop_objs);
         u=zeros(1,Global.M);
         deta=0;
        for i=1:Global.M
            for j=1:Global.N
                u(i)=u(i)+abs((Pop_objs(j,i)-Pre_Pop_objs(j,i))/(Max(i)-Min(i)));
            end
        end
        u=u/Global.N;
        for i=1:Global.M
            for j=1:Global.N
                f=(Pop_objs(j,i)-Pre_Pop_objs(j,i))/(Max(i)-Min(i));
                deta=deta+abs(f-u(i));
            end
        end
        I=deta/Global.M/Global.N;
end

