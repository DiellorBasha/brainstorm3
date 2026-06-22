function P = sa_crossing_profile(F, pairs, nProfiles, J)
% SA_CROSSING_PROFILE: Build sulcal-crossing profiles for Component 1's
% illustrative figure. For up to nProfiles facing-wall pairs, the surface mesh
% shortest path between the two walls descends through the fundus and back up;
% along it the normal turns ~pi while the Fiedler frame turns smoothly.
%
% USAGE:  P = sa_crossing_profile(F, pairs, nProfiles)
%         P = sa_crossing_profile(F, pairs, nProfiles, J)   % adds constrained scalar
%
%   F        : struct from sa_frames (uses Vtx, Nv, VertConn, fiedler).
%   pairs    : [nPairs x 2] facing-wall pairs (from sa_sulcal_walls).
%   nProfiles: max number of profiles to build.
%   J        : optional [nV x 3] current; if given, P(k).s is J.n along the path.
%
% OUTPUT P: struct array with fields path (vertex idx [L x 1]), arclen [L x 1],
% nAngle [L x 1] (cumulative normal turning, rad), fAngle [L x 1] (cumulative
% Fiedler-frame turning, rad), and s [L x 1] (constrained scalar, or []).

if nargin < 3 || isempty(nProfiles), nProfiles = 3; end
if nargin < 4, J = []; end
P = struct('path',{},'arclen',{},'nAngle',{},'fAngle',{},'s',{});
if isempty(pairs), return; end

Vtx = F.Vtx; Nv = F.Nv; fe1 = F.fiedler.e1;
% Weighted surface graph (Euclidean edge lengths).
[ei, ej] = find(triu(F.VertConn > 0, 1));
w = sqrt(sum((Vtx(ei,:) - Vtx(ej,:)).^2, 2));
G = graph(ei, ej, w, size(Vtx,1));

k = min(nProfiles, size(pairs,1));
for p = 1:k
    pth = shortestpath(G, pairs(p,1), pairs(p,2));
    if numel(pth) < 3, continue; end
    pth = pth(:);
    seg = sqrt(sum(diff(Vtx(pth,:),1,1).^2, 2));
    arclen = [0; cumsum(seg)];
    nAng = local_cum_turn(Nv(pth,:));
    fAng = local_cum_turn(fe1(pth,:));
    rec = struct('path', pth, 'arclen', arclen, 'nAngle', nAng, 'fAngle', fAng, 's', []);
    if ~isempty(J), rec.s = sum(J(pth,:) .* Nv(pth,:), 2); end
    P(end+1) = rec; %#ok<AGROW>
end
end

function ang = local_cum_turn(V)
% Cumulative absolute turning angle of a sequence of vectors (rad).
L = size(V,1);
ang = zeros(L,1);
for t = 2:L
    a = V(t-1,:); b = V(t,:);
    na = norm(a); nb = norm(b);
    if na < 1e-9 || nb < 1e-9, ang(t) = ang(t-1); continue; end
    c = max(min(dot(a,b)/(na*nb), 1), -1);
    ang(t) = ang(t-1) + acos(c);
end
end
