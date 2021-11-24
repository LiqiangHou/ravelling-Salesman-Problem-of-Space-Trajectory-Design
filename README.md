# ravelling-Salesman-Problem-of-Space-Trajectory-Design
 A conituous mapping technique is used to reformulate the mixed intger
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
% Author: 
%  Liqiang-Hou, School of Aeronautics and Astronautics, Shanghai Jiaotong University
% For any questions about the program, please email houliqiang@sjtu.edu.cn
% All right revserved
