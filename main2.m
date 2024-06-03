clc;
clear;
number_fun=1;
MaxIt=1000;
nPop=80;
[dim,low,high,Vio,GloMin,fun] = ProbInfo(number_fun);
for i=1:50
    disp(i);    
[value,fun_hist]=BES3_eng(nPop,MaxIt,low,high,dim,fun,Vio);
end
%end
%plot(fun_hist,'-','Linewidth',1.5)
%xlabel('Iteration')
%ylabel('fitness')
%legend('BES')
for j = 1:dim
        disp(['x', num2str(j), ': ', num2str(value.pos(j))]);
    end
