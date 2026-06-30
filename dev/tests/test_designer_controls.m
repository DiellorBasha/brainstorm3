% test_designer_controls - the 4 js kernels are discoverable, in the dynamic group, and realisable
ks = bst_eigfilter_kernel('list');
for n = {'gabor','travwave','resonator','stmatern'}
    assert(any(strcmp(ks, n{1})), sprintf('%s registered', n{1}));
    m = bst_eigfilter_kernel('info', n{1});
    assert(isfield(m,'domain') && strcmpi(m.domain,'js'), sprintf('%s domain js', n{1}));
end
% bandpass flags (gabor/travwave/resonator = wavelets; stmatern = low-pass)
assert(bst_eigfilter_kernel('info','gabor').bandpass==1, 'gabor bandpass');
assert(bst_eigfilter_kernel('info','stmatern').bandpass==0, 'stmatern lowpass');
disp('OK');
