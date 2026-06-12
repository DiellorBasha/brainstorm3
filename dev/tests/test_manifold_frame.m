function test_manifold_frame()
% TEST_MANIFOLD_FRAME  Headless regression for view_manifold's DeriveVertexFrame.
% Requires Brainstorm on path so view_manifold('DeriveVertexFrame', ...) dispatches.
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    % 2 hemispheres, nVh=2 each. Global vertices: L->[1;2], R->[3;4]; nVert=4.
    % Per vertex, grid = e1 + 1i*e2 for chosen orthonormal (e1,e2); rot = 1.
    %   L v1: e1=[1 0 0] e2=[0 1 0] -> N=[0 0 1]
    %   L v2: e1=[0 1 0] e2=[0 0 1] -> N=[1 0 0]
    %   R v3: e1=[0 0 1] e2=[1 0 0] -> N=[0 1 0]
    %   R v4: e1=[1 0 0] e2=[0 0 1] -> N=[0 -1 0]
    Emb(1).GlobalVertices = [1;2];
    Emb(1).vertex.grid     = [1 0 0; 0 1 0] + 1i*[0 1 0; 0 0 1];
    Emb(1).vertex.position = [10 0 0; 11 0 0];
    Emb(2).GlobalVertices = [3;4];
    Emb(2).vertex.grid     = [0 0 1; 1 0 0] + 1i*[1 0 0; 0 0 1];
    Emb(2).vertex.position = [0 10 0; 0 11 0];
    Ga(1).vertex.rotation = [1;1];
    Ga(1).singularity.vertices = uint32(1);          % local -> global 1
    Ga(2).vertex.rotation = [1;1];
    Ga(2).singularity.vertices = uint32([]);

    G = view_manifold('DeriveVertexFrame', Emb, Ga, 4);
    [nPass,nFail] = chk('U v1', isequal(G.U(1,:),[1 0 0]), nPass,nFail);
    [nPass,nFail] = chk('V v1', isequal(G.V(1,:),[0 1 0]), nPass,nFail);
    [nPass,nFail] = chk('N v1 = cross(U,V)', max(abs(G.N(1,:)-[0 0 1]))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('N v2', max(abs(G.N(2,:)-[1 0 0]))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('U v3 (hemi R scatter)', isequal(G.U(3,:),[0 0 1]), nPass,nFail);
    [nPass,nFail] = chk('anchor P v3 = position', isequal(G.P(3,:),[0 10 0]), nPass,nFail);
    [nPass,nFail] = chk('orthonormal U.V=0', max(abs(sum(G.U(1:4,:).*G.V(1:4,:),2)))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('unit |U|', max(abs(sqrt(sum(G.U(1:4,:).^2,2))-1))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('singularity globalized', isequal(G.Sing(:), 1), nPass,nFail);

    % off-support vertices stay zero (nVert larger than mapped indices)
    G2 = view_manifold('DeriveVertexFrame', Emb, Ga, 6);
    [nPass,nFail] = chk('off-support U zero', isequal(G2.U(5,:),[0 0 0]) && isequal(G2.U(6,:),[0 0 0]), nPass,nFail);

    % rotation rotates the in-plane frame: rot=1i -> U=-e2, V=e1 on L v1
    Ga2 = Ga; Ga2(1).vertex.rotation = [1i;1i];
    Gr = view_manifold('DeriveVertexFrame', Emb, Ga2, 4);
    [nPass,nFail] = chk('rot=1i -> U=-e2', max(abs(Gr.U(1,:)-[0 -1 0]))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('rot=1i -> V=e1',  max(abs(Gr.V(1,:)-[1 0 0]))<1e-12, nPass,nFail);

    % shape mismatch errors
    err=false; Bad=Emb; Bad(1).vertex.grid=[1 0 0]; try, view_manifold('DeriveVertexFrame', Bad, Ga, 4); catch, err=true; end
    [nPass,nFail] = chk('shape mismatch errors', err, nPass,nFail);

    fprintf('\n==== test_manifold_frame: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_manifold_frame: %d test(s) FAILED.', nFail); end
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
