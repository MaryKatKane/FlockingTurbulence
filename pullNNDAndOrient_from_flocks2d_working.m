%% Pull Nearest Neighbor Distances

% Goal: Pull nearest neighbor distances of particles, pull number of
% neighbors within 2x nearest neighbor distance and identify number of
% particles within, as well as average and minimum nearest neighbor
% distances in an attempt to find flocks, determine differences in theta
% between nearest neighboring particles.

% Reason for nearest neighbors - determine how close particles get and
% orient towards one another.

% Rename matrices for ease of coding:
Xp_vals = Xp_ov_V5taup001;
Yp_vals = Yp_ov_V5taup001;
Theta_vals = Theta_ov_V5taup001;

numTimePoints = length(Theta_vals);
numIndividuals = length(Theta_vals(:,1);

% Structure to save nearest neighbors
NearNeighStruct_V5taup001 = struct;

tempStruct = struct;

for ii = 1:numTimePoints % for each time point...
    
    for jj = 1:numIndividuals % for each individual...
        
        clearvars diffNeighDist diffNeighOrient minNeighID_win2x minNeighDist_win2x numNeigh_win2x minNeigh2xOrient
        
        % Pull X, Y, and Theta for individual at specific time point:
        tempOne_X = Xp_vals(jj,ii);
        tempOne_Y = Yp_vals(jj,ii);
        tempOne_Theta = Theta_vals(jj,ii);
        
        for kk = 1:length(Theta_vals(:,1)) % Checking individual against other individuals....
            if kk == jj % if it is the same time point...
                diffNeigh(kk) = NaN;
                diffOrient(kk) = NaN;
            elseif kk ~= jj % not the same particle...
                tempNeigh_X = Xp_vals(kk,ii);
                tempNeigh_Y = Yp_vals(kk,jj);
                tempNeigh_Theta = Theta_vals(kk,jj);
                
                diffNeigh(kk) = sqrt((tempNeigh_X - tempOne_X)^2 + (tempNeigh_Y - tempOne_Y)^2);
                diffOrient(kk) = tempNeigh_Theta - tempOne_Theta;
            end
            
            % find the minimum distance and corresponding difference in
            % Theta of the particles
            minNeighDist = min(diffNeigh);
            minNeighID = find(diffNeigh == minNeighDist);
            minNeighOrient = diffOrient(minNeighID);
            
            % find the number of neighbors within 2x the minimum neighbor
            % distance and the corresponding differences in orientation
            minNeighDist_2x = minNeighDist * 2;
            minNeighID_win2x = find(diffNeigh <= minNeighDist_2x);
            minNeighDist_win2x = diffNeigh(minNeighID_win2x);
            numNeigh_win2x = numel(find(diffNeigh <= minNeighDist_2x));
            minNeigh2xOrient = diffOrient(minNeighID_win2x);
            
            tempStruct(jj,ii).minNeighDist = minNeighDist;
            tempStruct(jj,ii).minNeighDiffOrient = minNeighOrient;
            tempStruct(jj,ii).minNeigh_wIn2x = minNeighDist_win2x;
            tempStruct(jj,ii).minNeighDiffOrient_wIn2x = minNeigh2xOrient;
            
        end
    end
end
            
NearNeighStruct_V5taup001 = tempStruct;            
            
            