function [BestSol, Convergence_curve ,timep]=BES3_eng(nPop,MaxIt,low,high,dim,fobj,VioFactor)
%nPop: size of population 
%MaxIt:number of iterations 
%low, high : space of Decision variables
%dim : number of Decision variables
%fobj : funcation 
 % paper citation : Alsattar, H. A., Zaidan, A. A., & Zaidan, B. B. (2020). Novel meta-heuristic bald eagle search optimisation algorithm. Artificial Intelligence Review, 53(3), 2237-2264.?
st=cputime;
% Initialize Best Solution
BestSol.cost = inf;
for i=1:nPop
   pop.pos(i,:) = low+(high-low).*rand(1,dim);
     [pop.cost(i),g,h]=fobj(pop.pos(i,:));
     v=sum(VioFactor.*max(0,[g,h]));
     pop.cost(i)=pop.cost(i)+v;
     pop.best_position(i,:) = pop.pos(i,:);
     pop.velocity(i,:) = zeros(1,dim);
    if pop.cost(i) < BestSol.cost
        BestSol.pos = pop.pos(i,:);
        BestSol.cost = pop.cost(i);
    end
end
 disp(num2str([0 BestSol.cost]))
for t=1:MaxIt
    %%               1- select_space 
    [pop BestSol s1(t)]=select_space(fobj,pop,nPop,BestSol,low,high,dim,VioFactor);
      %%                Update eagle velocity
    update_velocity(pop, BestSol, 0.7298, 1.4962, 2.0,nPop,dim);  % Example coefficients (c1, c2, learning rate)
    
    %%                Update eagle positions using velocity
    update_position(pop, low, high,nPop);
    %%                2- search in space
    [pop BestSol s2(t)]=search_space(fobj,pop,BestSol,nPop,low,high,VioFactor);
    %%                3- swoop
  [pop BestSol s3(t)]=swoop(fobj,pop,BestSol,nPop,low,high,dim,VioFactor);
  %testSum=testSum+BestSol.cost;     
  Convergence_curve(t)=BestSol.cost;
     disp(num2str([t BestSol.cost]));
    ed=cputime;
    timep=ed-st;
end
function [pop BestSol s1]=select_space(fobj,pop,npop,BestSol,low,high,dim,VioFactor)
Mean=mean(pop.pos);
% Empty Structure for Individuals
empty_individual.pos = [];
empty_individual.cost = [];
lm= 2;
s1=0;
for i=1:npop
    newsol=empty_individual;
    newsol.pos= BestSol.pos+ lm*rand(1,dim).*(Mean - pop.pos(i,:));
    newsol.pos = max(newsol.pos, low);
    newsol.pos = min(newsol.pos, high);
    [newsol.cost,g,h]=fobj(newsol.pos);
    v=sum(VioFactor.*max(0,[g,h]));
    newsol.cost=newsol.cost+v;
    if newsol.cost<pop.cost(i)
       pop.pos(i,:) = newsol.pos;
       pop.cost(i)= newsol.cost;
       s1=s1+1;
         if pop.cost(i) < BestSol.cost
          BestSol.pos= pop.pos(i,:);
          BestSol.cost=pop.cost(i); 
         end
    end
end
function update_velocity(pop, BestSol, w, c1, c2,npop,dim)
  % Update velocity based on social and cognitive components
  for i=1:npop
    social_term = c1 * rand(1,dim) .* (BestSol.pos - pop.pos(i,:));
    cognitive_term = c2 * rand(1,dim) .* (BestSol.pos - pop.pos(i,:));
    % Optional: Update learning rate over time (uncomment if using)
    % pop.velocity(i,:) = w * pop.velocity(i,:) + social_term + cognitive_term;
    pop.velocity(i,:) = social_term + cognitive_term;
  end


function update_position(pop, low, high,npop)
  % Update position using velocity
  for i=1:npop
    pop.pos(i,:) = pop.pos(i,:) + pop.velocity(i,:);
    pop.pos(i,:) = max(pop.pos(i,:), low);
    pop.pos(i,:) = min(pop.pos(i,:), high);
  end

function [pop best s1]=search_space(fobj,pop,best,npop,low,high,VioFactor)
Mean=mean(pop.pos);
a=10;
R=1.5;
% Empty Structure for Individuals
empty_individual.pos = [];
empty_individual.cost = [];
s1=0;
for i=1:npop-1
    A=randperm(npop);
pop.pos=pop.pos(A,:);
pop.cost=pop.cost(A);
        [x y]=polr(a,R,npop);
    newsol=empty_individual;
   Step = pop.pos(i,:) - pop.pos(i+1,:);
   Step1=pop.pos(i,:)-Mean;
   newsol.pos = pop.pos(i,:) +y(i)*Step+x(i)*Step1;
    newsol.pos = max(newsol.pos, low);
    newsol.pos = min(newsol.pos, high);
    [newsol.cost,g,h]=fobj(newsol.pos);
    v=sum(VioFactor.*max(0,[g,h]));
    newsol.cost=newsol.cost+v;
    if newsol.cost<pop.cost(i)
       pop.pos(i,:) = newsol.pos;
       pop.cost(i)= newsol.cost;
              s1=s1+1;
                  if pop.cost(i) < best.cost
                 best.pos= pop.pos(i,:);
                best.cost=pop.cost(i); 
            end
    end
end

function [pop best s1]=swoop(fobj,pop,best,npop,low,high,dim,VioFactor)
Mean=mean(pop.pos);
a=10;
R=1.5;
% Empty Structure for Individuals
empty_individual.pos = [];
empty_individual.cost = [];
s1=0;
for i=1:npop
     A=randperm(npop);
pop.pos=pop.pos(A,:);
pop.cost=pop.cost(A);
        [x y]=swoo_p(a,R,npop);
    newsol=empty_individual;
   Step = pop.pos(i,:) - 2*Mean;
   Step1= pop.pos(i,:)-2*best.pos;
   % Lévy Flight Parameters (adjust as needed)
   alpha = 1.5;  % Controls distribution shape
   beta = 1.5;   % Controls jump magnitude

   % Generate Lévy Jump using inverse transform method
   u = rand();
   v = rand();
   levy_jump = (sin(pi*alpha*(1 + beta)) * exp((alpha - 1)*log(v))) / (beta * power(u, (1/beta)));
   if u < 0.5
    levy_jump = -levy_jump;
   end

% Combine Step, Step1, and Lévy Jump for new position
newsol.pos = rand(1,length(Mean)).*best.pos + x(i)*Step + y(i)*Step1 + levy_jump;
% Random Walk step size (adjust as needed)
step_size = 0.1;  % Adjust step size as needed
random_walk = pop.pos(i, :) + step_size * (rand(1, dim) - 0.5);
% Combine Random Walk and Lévy Jump
newsol.pos = random_walk + levy_jump;    
% Evaluate new position with Lévy Jump
[newsol.cost,g,h]=fobj(newsol.pos);
    v=sum(VioFactor.*max(0,[g,h]));
    newsol.cost=newsol.cost+v;

    % Update eagle position if new position is better
    if newsol.cost < pop.cost(i)
        pop.pos(i,:) = newsol.pos;
        pop.cost(i) = newsol.cost;
             s1=s1+1;
                  if pop.cost(i) < best.cost
                 best.pos= pop.pos(i,:);
                best.cost=pop.cost(i); 
                  end
    end
end
function [xR yR]=swoo_p(a,R,N)
th = a*pi*exp(rand(N,1));
r  =th; %R*rand(N,1);
xR = r.*sinh(th);
yR = r.*cosh(th);
 xR=xR/max(abs(xR));
 yR=yR/max(abs(yR));
 
 function [xR yR]=polr(a,R,N)
%// Set parameters
th = a*pi*rand(N,1);
r  =th+R*rand(N,1);
xR = r.*sin(th);
yR = r.*cos(th);
 xR=xR/max(abs(xR));
 yR=yR/max(abs(yR));