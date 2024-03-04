% Files to use for optimizations at different temperatures
% Temp:  -25  -15  -05   05   15   25   35   45
% --------------------------------------------------------
% A123:   10   10   30   45   50   50   50   50
% ATL:    <2    5   15   20   30   40   50   50
% E1:      2    5   10   20   25   35   45   50
% E2:      2    5    5   15   25   35   45   50
% FORD53: 10   25   25   35   40   40   40
% FORD54: 10   25   25   35   40   40   40
% P14:     4   10   20   35   50   50   50   50
% SAM:     2    2    5   10   10   15   15   15
% 
% Note that the numeric indicator indicates the maximum C-rate used 
% in that test (e.g., "40" = 4.0C). Because of higher resistance at
% colder temperatures, we must scale down the maximum current to avoid
% exceeding voltage limits of the cell.
clear all
cellIDs = {'P14'};
temps = [  5   25   45];
mags =  {[30   50   50]}; % P14
         
