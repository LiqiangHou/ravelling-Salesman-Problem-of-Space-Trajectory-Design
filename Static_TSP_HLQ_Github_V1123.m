function Static_TSP_HLQ_Github_V1123
%) A conituous mapping technique is used to reformulate the mixed intger
% TSP into a continuous optimization probelm. A set of continuous
% paramteres, expected value and variance  of x and y to be traveled to the
% nex point are used to instead of the integer ID of point to be optimized.
% The program select the point to visit using priori expected value and
% varaiance, togther with the variance progagated from the previous point.
% Variance of the distance are then predicted using the updated variance at
% the current point and the previous point. After cntinuous reformulation,
% a gradient-based optimizer is used to search the optimal solution. Three
% bencmarks of TSPLIB, burma14,bays29 and korA100 are used to test the
% algorithm.
%

% Benchmarks are from:http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/ 

clear all;
close all;
format short g;
global graph

testcase = 1;

% parameter  setting

if(testcase == 1)
    graph    = construct_graph_burma14();
    ID0      = 13;
    sigma_x0 = [0.5, 0.5];
    
    
    init_step = [ graph.min_abs_dst_x /(graph.max_dx - graph.min_dx) , ...
                  graph.min_abs_dst_y /(graph.max_dy - graph.min_dy) ] ;
    mu_x0     = 0.5 + init_step;
elseif(testcase == 2)
    graph    = construct_graph_bay29();
    ID0      = 1;
    sigma_x0 = [1.0, 1.0];   

    
    init_step = [ graph.min_abs_dst_x /(graph.max_dx - graph.min_dx) , ...
                  graph.min_abs_dst_y /(graph.max_dy - graph.min_dy) ] ;
    mu_x0     = 0.5 + init_step;
else
    graph    = construct_graph_kroA100();
    ID0      = 1;

    
    sigma_x0 = [1.0, 1.0];    
    
    init_step = [ graph.min_abs_dst_x /(graph.max_dx - graph.min_dx) , ...
                  graph.min_abs_dst_y /(graph.max_dy - graph.min_dy) ] ;
    mu_x0     = 0.5 - init_step;
end
 
 
% lower bound and upper bound 
[lb,ub] = Initialize_design_varaibles(graph);


%
x0    = [];
for i=1:graph.n

    x0 = [x0,   ...
          mu_x0 ,...            
          sigma_x0];
end
% 

% SQP is used to search the optimal route 
fun     = @(x)my_route_cost(x,graph,ID0);
options = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',1.0e-30,'StepTolerance',1.0e-10,'TolFun',1.0e-30,'TolX',1.0e-15,'MaxIter',800,'Display','iter','OutputFcn',@outfun);

[fval] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
 

return
end



function J = my_route_cost(x,graph,ID0)
global tour_fitness_sum_Penalty_Log_P

    
sum_Penalty_Log_P = 0;              
 
 
 
% remove ID0 from the list
id_pre       = ID0;
node_list    = graph.ID;
[node_list]  = remove_visited_node(ID0,node_list);

SIGMA_x_pre = zeros(2,2);
for i_node = 1:graph.n
    

%         
    [mu_x,SIGMA_x] = design_varaibles_of_node(x,i_node,SIGMA_x_pre,graph);
    
    % estimated y
    mu_y = sqrt(mu_x(1)^2 + mu_x(2)^2);
     if(i_node == graph.n)
         node_list = ID0;
     end
        % likelihood of candidate
        [id_likelihood_xlist_ylist, SIGMA_xy, SIGMA_yx, SIGMA_y]  = Probilities_of_Node(graph,id_pre,node_list,mu_x,SIGMA_x);
        
        % select
        [id,logLikelihood,post_SIGMA_x] = Select_Next_Node(id_likelihood_xlist_ylist,mu_x,mu_y,SIGMA_x, SIGMA_xy, SIGMA_yx, SIGMA_y);
        
        selected_id(i_node) = id;
        y_hat_set(i_node)   =  mu_y;                        


        
         % update list of candidtae nodes 
        [node_list] = remove_visited_node(id,node_list);
        sum_Penalty_Log_P     = sum_Penalty_Log_P + logLikelihood; 
        
        SIGMA_x_pre = post_SIGMA_x;
        id_pre      = id;
        

end


% fiteness function
whole_route      = [ID0,selected_id];
fitness = fitnessFunction ( whole_route , graph);
 
sum_y_hat    = sum(y_hat_set);

%           ;
J         = sum_y_hat...
          + sum_Penalty_Log_P...
          ;     
tour_fitness_sum_Penalty_Log_P = [whole_route,fitness,J];


return
end
%
function [F] = mean_and_cov_of_individual_node(x)
%
dJ_dx = x(1)/norm(x);
dJ_dy = x(2)/norm(x);
 
 
% first order Jacobian of cost function (eulidian distance is used)
F       = [dJ_dx,dJ_dy];


return
end
 
%

 
%
%
function [id_likelihood_xlist_ylist, SIGMA_xy, SIGMA_yx, SIGMA_y] = Probilities_of_Node(graph,id_pre,node_list,mu_x,SIGMA_x)
global V %  
 
%
n_node     = length(node_list);

V = (graph.min_edge )^2;                                   


mu_node_pre = [graph.node(id_pre).x , graph.node(id_pre).y] ;



% predict value and cov of cost. F: Jacobian matrix of individual cost of the node
[F]      = mean_and_cov_of_individual_node(mu_x)  ;
SIGMA_xy = SIGMA_x*F';
SIGMA_yx = F*SIGMA_x;
SIGMA_y  = F*SIGMA_x*F' + V;

for i =1:n_node
    x  = [graph.node(node_list(i)).x , graph.node(node_list(i)).y] - mu_node_pre;
    y  = sqrt(x(1)^2 + x(2)^2);

    
    x_list(i,:) = x;
    y_list(i)   = y;

   
    logLikelihood(i) = [(x - mu_x)]*invChol_mex(SIGMA_x )*[(x - mu_x) ]' + log(det(SIGMA_x));


end

id_likelihood_xlist_ylist = [node_list',logLikelihood',x_list,y_list'];
 

return
end
 


 
function [id,logLikelihood,post_SIGMA_x]= Select_Next_Node(id_likelihood_xlist_ylist,mu_x,mu_y,SIGMA_x, SIGMA_xy, SIGMA_yx, SIGMA_y)
global V %  
[Max_L,I]   = min(id_likelihood_xlist_ylist(:,2));
%
id      = id_likelihood_xlist_ylist(I,1);
y       = id_likelihood_xlist_ylist(I,5);

% posterior covariance 
post_SIGMA_x  =  SIGMA_x - SIGMA_xy * (SIGMA_y + V)^(-1) * SIGMA_yx;

% log-liklihood of cost function
logLikelihood = (y - mu_y) * (SIGMA_y + V)^(-1) * (y - mu_y)' + log(det(SIGMA_y + V));

 

return
end
%
%
%
function [mu_x,SIGMA_x,kappa] = design_varaibles_of_node(x,i_node,SIGMA_x_pre,graph)
global W 
    
    n_var        = 4;
    
    m_target     = x((i_node - 1)*n_var + 1  : (i_node - 1)*n_var + 2);         % dx,dy
    d_target     = x((i_node - 1)*n_var + 3  : (i_node - 1)*n_var + 4);         % std_x,std_y
    


 

    min_d_pos     = [graph.min_dx        graph.min_dy];
    max_d_pos     = [graph.max_dx        graph.max_dy];
    max_std_dr    = [graph.max_abs_dst_x  ,...
                     graph.max_abs_dst_y  ];      

    %
    mu_x     = min_d_pos + (max_d_pos - min_d_pos).*m_target;   % from [0  1] to [min_dr max_dr]
    
    W = [   (graph.min_abs_dst_x)^2           0
             0                               (graph.min_abs_dst_y)^2 ]; % 9425.7
    SIGMA_x  = diag(d_target .* max_std_dr ) + W + SIGMA_x_pre;


return
end
 
function [new_node_list] = remove_visited_node(ID,node_list)
 
new_node_list = node_list;
 
% remove the node that has been visted
IID = find(node_list == ID);
new_node_list(IID) = [];
 
 
 
return
end
 
 
 
 
function [lb,ub] = Initialize_design_varaibles(graph)
% the number of edges is n-1 
lb = [];
ub = [];
 
 
    
for i_node = 1:graph.n
    
    [lb_node,ub_node] = Initialize_design_varaibles_node();
    
    lb = [lb,lb_node];
    ub = [ub,ub_node];
end
 
 
 
return
end
 
 
function [lb,ub] = Initialize_design_varaibles_node()

%

min_dx = 0.0;
min_dy = 0.0;
 
max_dx = 1.0;
max_dy = 1.0;
 

min_std_x = 1.0e-6;
min_std_y = 1.0e-6;

max_std_x = 1.0;
max_std_y = 1.0;
 

 
 
% expected dx and dy 
lb_dx =  min_dx;
ub_dx =  max_dx;
 
lb_dy =  min_dy;
ub_dy =  max_dy;
 
% std of dx and dy
lb_std_x = min_std_x;
ub_std_x = max_std_x;
 
lb_std_y = min_std_y;
ub_std_y = max_std_y;

 
%

lb(1) = lb_dx;
lb(2) = lb_dy;
lb(3) = lb_std_x;
lb(4) = lb_std_y;

%
ub(1) = ub_dx;
ub(2) = ub_dy;
ub(3) = ub_std_x;
ub(4) = ub_std_y;

 
return
end
%---------------------------------------------
 
 
 
function [ fitness ] = fitnessFunction ( tour , graph)
 
 
fitness = 0;
 
for i = 1 : length(tour) -1
    
    currentNode = tour(i);
    nextNode = tour(i+1);
    
    fitness = fitness + graph.edges( currentNode ,  nextNode );
    
end
return
end
 
 
 
 
 
 
 
 
function graph = construct_graph_bay29()
% bay29 nodes

DISPLAY_DATA_SECTION = [
   1    1150.0  1760.0
   2     630.0  1660.0
   3      40.0  2090.0
   4     750.0  1100.0
   5     750.0  2030.0
   6    1030.0  2070.0
   7    1650.0   650.0
   8    1490.0  1630.0
   9     790.0  2260.0
  10     710.0  1310.0
  11     840.0   550.0
  12    1170.0  2300.0
  13     970.0  1340.0
  14     510.0   700.0
  15     750.0   900.0
  16    1280.0  1200.0
  17     230.0   590.0
  18     460.0   860.0
  19    1040.0   950.0
  20     590.0  1390.0
  21     830.0  1770.0
  22     490.0   500.0
  23    1840.0  1240.0
  24    1260.0  1500.0
  25    1280.0   790.0
  26     490.0  2130.0
  27    1460.0  1420.0
  28    1260.0  1910.0
  29     360.0  1980.0];

% ID of the nodes
graph.ID = DISPLAY_DATA_SECTION(:,1)';

% pos of nodes
x = DISPLAY_DATA_SECTION(:,2);
y = DISPLAY_DATA_SECTION(:,3);
 
graph.center.x = x(1);
graph.center.y = y(1);
%
graph.n = length(x);


for i = 1 : graph.n
    graph.node(i).x = x(i);
    graph.node(i).y = y(i);
end
 
graph.edges  = zeros(  graph.n , graph.n );
graph.dist_x = zeros(  graph.n , graph.n );
graph.dist_y = zeros(  graph.n , graph.n );
 


for i = 1 : graph.n
    for j = 1: graph.n
        x1 = graph.node(i).x ;
        x2 = graph.node(j).x;
        y1 = graph.node(i).y;
        y2 = graph.node(j).y;
        
        graph.edges(i,j) = sqrt(  (x1 - x2) ^2 + (y1 - y2)^2  );
        graph.dist_x(i,j)  = (x1 - x2);
        graph.dist_y(i,j)  = (y1 - y2);
        
    end
end

%
edges  = graph.edges;
dist_x = graph.dist_x;
dist_y = graph.dist_y;


graph.min_dx   = min(min(dist_x)) ;
graph.max_dx   = max(max(dist_x)) ;

graph.min_dy   = min(min(dist_y)) ;
graph.max_dy   = max(max(dist_y)) ;


abs_dist_x = abs(dist_x);
abs_dist_y = abs(dist_y);

graph.min_abs_dst_x = min(abs_dist_x(abs_dist_x > 0));
graph.min_abs_dst_y = min(abs_dist_y(abs_dist_y > 0));

graph.max_abs_dst_x = max(abs_dist_x(abs_dist_x > 0));
graph.max_abs_dst_y = max(abs_dist_y(abs_dist_y > 0));




graph.min_edge = min(edges(edges >0)) ;
graph.max_edge = max(edges(edges >0)) ;



graph.min_std_dist_x = min(std(dist_x)); 
graph.max_std_dist_x = max(std(dist_x)); 

graph.min_std_dist_y = min(std(dist_y)); 
graph.max_std_dist_y = max(std(dist_y)); 

graph.min_std_edges =  min(std(edges));
graph.max_std_edges =  max(std(edges));

return
end
 
 
 
function [ ] = drawBestTour(currentSolution , graph, fitness)
figure(2);
for i = 1 : length(currentSolution) - 1
    
    currentNode = currentSolution(i);
    nextNode =  currentSolution(i+1);
    
    x1 = graph.node(currentNode).x;
    y1 = graph.node(currentNode).y;
    
    x2 = graph.node(nextNode).x;
    y2 = graph.node(nextNode).y;
    
    X = [x1 , x2];
    Y = [y1, y2];
    
    plot (X, Y, '-r');
    text(x1+0.2, y1,num2str(currentSolution(i)));
    hold on;

end


title(['Best tour','Total length: ',num2str(fitness)]);
box('on');
hold off;

return
end
 
 



function stop = outfun(x,optimValues,state)
global tour_fitness_sum_Penalty_Log_P
global graph
 
stop = false;
 
switch state
    case 'init'
        hold on
        
    case 'iter'
        % Concatenate current point and objective function
        % value with history. x must be a row vector.
        hold off
 
        data_tour_fitness_sum_Penalty_Log_P = tour_fitness_sum_Penalty_Log_P;
        currentSolution = tour_fitness_sum_Penalty_Log_P(1:end-2);
        fitness         = tour_fitness_sum_Penalty_Log_P(end-1);
        drawBestTour(currentSolution , graph, fitness);
 
        outmsg = [ ' Shortest length = ' , num2str(fitness) ];
        disp(outmsg)
    case 'done'
        hold off
    otherwise
end
 
return
end
 
 


 
 
function graph = construct_graph_burma14()
% 14 nodes
x = [16.47000,16.47000,20.09000,22.39000,25.23000,22.00000,20.47000,17.20000,16.30000,14.05000,16.53000,21.52000,19.41000,20.09000];
y = [96.10000,94.44000,92.54000,93.37000,97.24000,96.05000,97.02000,96.29000,97.38000,98.12000,97.38000,95.59000,97.13000,94.55000];


% pos of vertex
graph.ID = [1:14];
 
%
graph.n = length(x);
 
for i = 1 : graph.n
    graph.node(i).x = x(i);
    graph.node(i).y = y(i);
end
 
graph.edges  = zeros(  graph.n , graph.n );
graph.dist_x = zeros(  graph.n , graph.n );
graph.dist_y = zeros(  graph.n , graph.n );
 


for i = 1 : graph.n
    for j = 1: graph.n
        x1 = graph.node(i).x ;
        x2 = graph.node(j).x;
        y1 = graph.node(i).y;
        y2 = graph.node(j).y;
        
        graph.edges(i,j) = sqrt(  (x1 - x2) ^2 + (y1 - y2)^2  );
        graph.dist_x(i,j)  = [(x1 - x2)];
        graph.dist_y(i,j)  = [(y1 - y2)];
        
    end
end
%
edges  = graph.edges;
dist_x = graph.dist_x;
dist_y = graph.dist_y;


graph.min_dx   = min(min(dist_x)) ;
graph.max_dx   = max(max(dist_x)) ;

graph.min_dy   = min(min(dist_y)) ;
graph.max_dy   = max(max(dist_y)) ;


abs_dist_x = abs(dist_x);
abs_dist_y = abs(dist_y);

graph.min_abs_dst_x = min(abs_dist_x(abs_dist_x > 0));
graph.min_abs_dst_y = min(abs_dist_y(abs_dist_y > 0));

graph.max_abs_dst_x = max(abs_dist_x(abs_dist_x > 0));
graph.max_abs_dst_y = max(abs_dist_y(abs_dist_y > 0));




graph.min_edge = min(edges(edges >0)) ;
graph.max_edge = max(edges(edges >0)) ;



graph.min_std_dist_x = min(std(dist_x)); 
graph.max_std_dist_x = max(std(dist_x)); 

graph.min_std_dist_y = min(std(dist_y)); 
graph.max_std_dist_y = max(std(dist_y)); 

graph.min_std_edges =  min(std(edges));
graph.max_std_edges =  max(std(edges));
return
end
 
 


 
 
function graph = construct_graph_kroA100()
% 100 nodes
nodes = [1380,2848,3510,457,3888,984,2721,1286,2716,738,1251,2728,3815,3683,1247,123,1234,252,611,2576,928,53,1807,274,2574,178,2678,1795,3384,3520,1256,1424,3913,3085,2573,463,3875,298,3479,2542,3955,1323,3447,2936,1621,3373,1393,3874,938,3022,2482,3854,376,2519,2945,953,2628,2097,890,2139,2421,2290,1115,2588,327,241,1917,2991,2573,19,3911,872,2863,929,839,3893,2178,3822,378,1178,2599,3416,2961,611,3113,2597,2586,161,1429,742,1625,1187,1787,22,3640,3756,776,1724,198,3950;939,96,1671,334,666,965,1482,525,1432,1325,1832,1698,169,1533,1945,862,1946,1240,673,1676,1700,857,1711,1420,946,24,1825,962,1498,1079,61,1728,192,1528,1969,1670,598,1513,821,236,1743,280,1830,337,1830,1646,1368,1318,955,474,1183,923,825,135,1622,268,1479,981,1846,1806,1007,1810,1052,302,265,341,687,792,599,674,1673,1559,558,1766,620,102,1619,899,1048,100,901,143,1605,1384,885,1830,1286,906,134,1025,1651,706,1009,987,43,882,392,1642,1810,1558];
x = nodes(1,:);
y = nodes(2,:);


%
graph.n = length(x);


% pos of vertex
graph.ID = [1:graph.n];
 
 
for i = 1 : graph.n
    graph.node(i).x = x(i);
    graph.node(i).y = y(i);
end
 
graph.edges  = zeros(  graph.n , graph.n );
graph.dist_x = zeros(  graph.n , graph.n );
graph.dist_y = zeros(  graph.n , graph.n );
 


for i = 1 : graph.n
    for j = 1: graph.n
        x1 = graph.node(i).x ;
        x2 = graph.node(j).x;
        y1 = graph.node(i).y;
        y2 = graph.node(j).y;
        
        graph.edges(i,j) = sqrt(  (x1 - x2) ^2 + (y1 - y2)^2  );
        graph.dist_x(i,j)  = [(x1 - x2)];
        graph.dist_y(i,j)  = [(y1 - y2)];
        
    end
end
%
edges  = graph.edges;
dist_x = graph.dist_x;
dist_y = graph.dist_y;


graph.min_dx   = min(min(dist_x)) ;
graph.max_dx   = max(max(dist_x)) ;

graph.min_dy   = min(min(dist_y)) ;
graph.max_dy   = max(max(dist_y)) ;


abs_dist_x = abs(dist_x);
abs_dist_y = abs(dist_y);

graph.min_abs_dst_x = min(abs_dist_x(abs_dist_x > 0));
graph.min_abs_dst_y = min(abs_dist_y(abs_dist_y > 0));

graph.max_abs_dst_x = max(abs_dist_x(abs_dist_x > 0));
graph.max_abs_dst_y = max(abs_dist_y(abs_dist_y > 0));




graph.min_edge = min(edges(edges >0)) ;
graph.max_edge = max(edges(edges >0)) ;



graph.min_std_dist_x = min(std(dist_x)); 
graph.max_std_dist_x = max(std(dist_x)); 

graph.min_std_dist_y = min(std(dist_y)); 
graph.max_std_dist_y = max(std(dist_y)); 

graph.min_std_edges =  min(std(edges));
graph.max_std_edges =  max(std(edges));
return
end
 
 

 