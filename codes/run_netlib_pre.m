%% netlib
Probname = {
'sebapre.mat','sc50bpre.mat','sc50apre.mat','afiropre.mat',...                                   % less than 1KB  4
'sc105pre.mat','beaconfdpre.mat', 'vtp_basepre.mat','sc205pre.mat','recipepre.mat',...
'scagr7pre.mat','kb2pre.mat','scorpionpre.mat','stocfor1pre.mat', 'adlittlepre.mat',...
'share2bpre.mat','lotfipre.mat','standatapre.mat', 'standgubpre.mat', 'sctap1pre.mat',...
'bore3dpre.mat','scagr25pre.mat','boeing2pre.mat','scsd1pre.mat','share1bpre.mat',...            % less than 5KB  20
'israelpre.mat','aggpre.mat','standmpspre.mat', 'brandypre.mat', 'scsd6pre.mat',...
'etamacropre.mat', 'gangespre.mat', 'degen2pre.mat', 'finnispre.mat','capripre.mat',...
'shellpre.mat','sctap2pre.mat','agg2pre.mat','agg3pre.mat',...                                   % less than 10KB 14
'bandmpre.mat','e226pre.mat','grow7pre.mat','scfxm1pre.mat','ship04spre.mat',...
'qap08pre.mat','sctap3pre.mat','boeing1pre.mat','scsd8pre.mat','fffff800pre.mat',...
'ship08spre.mat','scrs8pre.mat','stocfor2pre.mat', 'ship04lpre.mat','ship12spre.mat',...
'tuffpre.mat','bnl1pre.mat','pds-02pre.mat','modszk1pre.mat',...                                 % less than 20KB 19
'grow15pre.mat', 'marospre.mat','scfxm2pre.mat','fit1ppre.mat','czprobpre.mat',...
'ken-07pre.mat','stairpre.mat','ship08lpre.mat','fit1dpre.mat', 'grow22pre.mat',...
'scfxm3pre.mat',...                                                                              % less than 30KB 11
'woodwpre.mat','degen3pre.mat','peroldpre.mat','cre-cpre.mat','ship12lpre.mat',...
'pilot4pre.mat','bnl2pre.mat','cre-apre.mat',...                                                 % less than 40KB 8
'nesmpre.mat','25fv47pre.mat','pilot.wepre.mat','d6cubepre.mat','pilotnovpre.mat',...
'qap12pre.mat','pilot.japre.mat','cyclepre.mat','pds-06pre.mat',...                              % less than 100KB 9
'fit2ppre.mat','80bau3bpre.mat','wood1ppre.mat','greenbebpre.mat','greenbeapre.mat',...
'ken-11pre.mat', 'fit2dpre.mat','d2q06cpre.mat',...                                              % less than 150KB 8
'cre-dpre.mat','maros-r7pre.mat','qap15pre.mat','pds-10pre.mat','pilotpre.mat',...
'cre-bpre.mat', 'ken-13pre.mat','osa-07pre.mat','pilot87pre.mat','pds-20pre.mat',...             % less than 500KB  10
'osa-14pre.mat','ken-18pre.mat',...                                                              % less than 1MB   2
'osa-30pre.mat','osa-60pre.mat'                                                                  % more than 1MB   2
};   % total 107

%%
nlen = length(Probname);
problist = 1 : 24;
for dprob = 1:length(problist)
    pid  = problist(dprob);
    name = Probname{pid};
    load(strcat(filesep,name),'Model');
    out = preprocess(Model);
    c = out.c;
    A = out.A; 
    b = out.b;
    [m, n] = size(A);
    x0 = abs(randn(n, 1));
    disp(name);
    run Test_lp_problems.m
end





    
