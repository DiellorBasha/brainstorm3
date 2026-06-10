function test_dirac_eigenmode_field
% Pure reconstruction helper: J(c, vertex) = imag( Phi_D * gainRow(modes_h)' ).
% Synthetic, deterministic, no GUI / no nxr.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);

nCh = 5; vH1 = [1 2 3]; vH2 = [4 5]; nVert = 5; Kh = 2;
hemis = {vH1(:), vH2(:)};

DE = repmat(struct('Vectors',[],'GlobalVertices',[]), 1, 2);
for hh = 1:2
    nVh = numel(hemis{hh});
    raw = reshape(sin((1:4*nVh)' * (1:(Kh+2)) * pi/13 + hh), 4*nVh, Kh+2);  % full rank
    Phi = orth(raw); Phi = Phi(:, 1:Kh);
    DE(hh).Vectors = Phi; DE(hh).GlobalVertices = hemis{hh};
end
CompHM = struct('nModes', 2*Kh, 'ModeHemisphere', [1;1;2;2], ...
                'HemiGlobalVertices', {{vH1(:), vH2(:)}});

GainRows = reshape(cos((1:nCh*2*Kh) * 0.37), nCh, 2*Kh);   % deterministic coefficients

J = bst_dirac_eigenmode_field(DE, GainRows, CompHM);
assert(isequal(size(J), [nCh, 3*nVert]), 'J must be [m x 3*nVert].');

% Independent oracle: per channel, per hemisphere, imag-part of Phi*g_h at each vertex.
for c = 1:nCh
    for hh = 1:2
        vH = hemis{hh};
        cols = (CompHM.ModeHemisphere == hh);
        psi  = DE(hh).Vectors * GainRows(c, cols).';        % [4Vh x 1]
        for k = 1:numel(vH)
            v   = vH(k);
            ref = [psi(4*(k-1)+2), psi(4*(k-1)+3), psi(4*(k-1)+4)];
            got = J(c, 3*(v-1)+(1:3));
            assert(max(abs(got - ref)) < 1e-10, sprintf('recon mismatch ch %d vert %d', c, v));
        end
    end
end

% Linearity sanity: scaling a channel's coefficients scales its field.
J2 = bst_dirac_eigenmode_field(DE, 2*GainRows, CompHM);
assert(max(abs(J2(:) - 2*J(:))) < 1e-9, 'reconstruction must be linear.');

disp('ALL TESTS PASSED');
end
