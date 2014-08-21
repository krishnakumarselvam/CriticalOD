%List of Inputs
%1) Randonly chosen trial points
%2) TopODIndices
%3) baseODMatrix

%4) Queing information


clc;clear; 
load ('Inputs/RandomTrialPoints', 'randomTrialPoints');
load ('Inputs/RandomTrialPoints', 'TopODIndices');
baseODMatrix = textread('Inputs/ODpairs.txt'); 




























%%
% function [c, ceq] = nonlcon(x0)
%        %Non linear inequalities
%         c=[];
%        %Non linear equalities
%            % 1) Relating the OD matrix to the arrival rates
%            % 2) The queuing equations
%    
%         % 1) Relating the OD matrix to the arrival rates
%         %---------------------------------------------------
%         
%         %First get the full OD Vector as a function of x0. The OD vector is
%         %the third column in the baseODMatrix variable (Origin,
%         %Destination, Demand). That is, it is a vector of the demand
%         %quantities
%         ODVector = baseODMatrix(:,3);
%         ODVector (TopODIndices) = x0;
%         
%         %Now connect the arrival rates in different queues to the OD matrix
%         %Most Arrival rates are constant - some would depend on the values
%         %in x0. 
%         
%         % Gamma (of a few indices) = Linear function of x0
%         
%         % 2) The queuing equations
%         %---------------------------------------------------
%       
%         % Write down the 3 variable queuing model equations between
%         % the queuing variables. 
%         
% 
%         
% 
% 
% end

