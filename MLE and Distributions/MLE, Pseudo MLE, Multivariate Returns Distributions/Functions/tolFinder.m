% Finding the optimal tolerance 
function tol = tolFinder(method, data, initvec)
% initialise
left = 1e-8;
right = 1e-1;
maxit=100;

% titer
if nargin < 3
    % accurate value
    [dfSol, muSol, cSol] = method(data, left); 
    sol = [dfSol, muSol, cSol];

    for i = 1:maxit
        % get middle tolerance
        mid = (left + right) / 2;
        
        % evaluate result at left and middle tolerance
        [dfL, muL, cL] = method(data, left);
        evalL = [dfL, muL, cL];
        [dfM, muM, cM] = method(data, mid);
        
        % if middle point is of four digit accuracy, set it as the new left
        % point, if it is not, set it as the new right point
        if (floor_sf(dfL, 4) == floor_sf(dfM, 4)) && (floor_sf(muL, 4) == ...
                floor_sf(muM, 4)) && (floor_sf(cL, 4) == floor_sf(cM, 4))
            left = mid;
        else
            right = mid;
        end
    end
    disp(sol);
else
    % accurate value
    sol = method(data, left, initvec)'; 

    for i = 1:maxit
        % get middle tolerance
        mid = (left + right) / 2;
        
        % evaluate result at left and middle tolerances
        evalL = method(data, left, initvec)';
        [dfL, muL, cL] = deal(evalL(1), evalL(2), evalL(3));
        evalM = method(data, mid, initvec)';
        [dfM, muM, cM] = deal(evalM(1), evalM(2), evalM(3));

        % if middle point is of four digit accuracy, set it as the new left
        % point, if it is not, set it as the new right point
        if (floor_sf(dfL, 4) == floor_sf(dfM, 4)) && (floor_sf(muL, 4) == ...
                floor_sf(muM, 4)) && (floor_sf(cL, 4) == floor_sf(cM, 4))
            left = mid;
        else
            right = mid;
        end
    end
    disp(sol);
end
disp(evalL); % TEST OUTPUT OPT TOL
tol = left;
disp(left);
end