function D = bst_time_derivative(F, dt, order)
% BST_TIME_DERIVATIVE  Backward finite-difference time derivative of a field.
%   F     : [n x m] consecutive frames, columns oldest -> newest (last col = current).
%           Must have at least order+1 columns; the newest order+1 are used.
%   dt    : time step (s)
%   order : 0 (field) | 1 (velocity) | 2 (acceleration)
%   D     : [n x 1]
% Author: Diellor Basha, 2026
    if size(F,2) < order+1
        error('bst_time_derivative:tooFewFrames', 'need %d frames for order %d (got %d)', order+1, order, size(F,2));
    end
    switch order
        case 0, D = F(:,end);
        case 1, D = (F(:,end) - F(:,end-1)) / dt;
        case 2, D = (F(:,end) - 2*F(:,end-1) + F(:,end-2)) / dt^2;
        otherwise, error('bst_time_derivative:badOrder', 'order must be 0, 1, or 2');
    end
end
