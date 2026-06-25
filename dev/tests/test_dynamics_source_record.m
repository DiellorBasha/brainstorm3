function test_dynamics_source_record
    % Prime a hidden figure with a fake DynamicsOverlay cache; OnRecord must read Ht.Curl.
    hF = figure('Visible','off');  c = onCleanup(@() close(hF));
    nV = 4; Ht = struct('Div',zeros(nV,1),'Curl',[0;5;-3;0],'Phi',zeros(nV,1),'Psi',zeros(nV,1));
    D = struct('Cov',[],'Op','Curl','Cache',containers.Map('KeyType','double','ValueType','any'), ...
               'srcDS',1,'srcResult',1,'iTess',1,'nV',nV);
    D.Cache(1) = Ht;  setappdata(hF,'DynamicsOverlay',D);
    s = view_dynamics('PickScalar', D.Cache(1), 'Curl');
    assert(isequal(s, Ht.Curl), 'overlay cache wiring: PickScalar(Curl) must read Ht.Curl');
    fprintf('PASS test_dynamics_source_record (cache+pick wiring)\n');
end
