function model_fit(subject_n)
%% function for model fit



% Repository paths
scriptFolder = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptFolder));

dataFolder     = fullfile(repoRoot,'data');
behaviorFolder = fullfile(dataFolder,'behavior');
modelFolder    = fullfile(dataFolder,'model');
resultsFolder = fullfile(modelFolder,'results');

addpath(genpath(fullfile(repoRoot,'functions','bads-master')));





for s=subject_n


    % subject s

    tic;

    % load data
    load(fullfile(behaviorFolder,'data.mat'))

    this_subject=[];
    % aggregate drug conditions
    for dd=1:3

        this_subject_id=cell2mat(All_Data(s,dd));
        this_subject=[this_subject;this_subject_id];

    end


    % extract relevant trial info

    % choices
    Resp_ses=this_subject(:,1);
    % contrasts
    All_contrasts=[this_subject(:,12:21),this_subject(:,2:11)];
    %RT
    All_RT=this_subject(:,24);
    %confiidence report
    All_conf=this_subject(:,22)+1;
    %generative mean
    All_gm=this_subject(:,23);
    % session index
    All_ses_id=this_subject(:,25:26);

    % take out nan vaues
    id_nan=(isnan(All_conf));
    Resp_ses(id_nan)=[];
    All_contrasts(id_nan,:)=[];
    All_RT(id_nan,:)=[];
    All_conf(id_nan,:)=[];
    All_gm(id_nan,:)=[];
    All_ses_id(id_nan,:)=[];





    global dat;
    dat.t=10;

    % choices
    dat.response=Resp_ses';

    % left right contrasts
    dat.right=All_contrasts(:,1:10)';
    dat.left=All_contrasts(:,11:20)';
    dat.rt=All_RT;
    dat.ses_id=All_ses_id;
    dat.conf=All_conf;



    % high/low confidence trials
    id_hc=find(All_conf==2);
    id_lc=find(All_conf==1);


    % CV mantaining the proportion oif High/Low conf trials across folds
    id_cv_lc=cvpartition(size(id_lc,1),'Kfold',5);
    id_cv_hc=cvpartition(size(id_hc,1),'Kfold',5);

    All_stimulus=[ ];
    All_left=[ ];
    All_right=[ ];
    All_ch_model=[ ];
    All_Conf=[];
    Real_data=[];
    LL_cv     = nan(5,1);
    All_param = nan(5,11);
    % 5 fold CV
    for cv=1:5

        tic
        cv

        % tak train and test Id
        id_train_hc=training(id_cv_hc,cv);
        id_test_hc=test(id_cv_hc,cv);
        id_train_lc=training(id_cv_lc,cv);
        id_test_lc=test(id_cv_lc,cv);



        id_train=[id_lc(id_train_lc);id_hc(id_train_hc)];
        id_test=[id_lc(id_test_lc);id_hc(id_test_hc)];



        % train trials
        global dat;
        dat.t=10;

        % choices
        dat.response=Resp_ses(id_train)';

        % left right contrasts
        dat.right=All_contrasts(id_train,1:10)';
        dat.left=All_contrasts(id_train,11:20)';
        dat.rt=All_RT(id_train);

        dat.conf=All_conf(id_train);
        dat.ses_id=All_ses_id(id_train,:);

        dat.gm=All_gm(id_train);

        % number of initial values
        n_iteration=500;


        leak_range=((1.9).*rand(2000,1)+0.04)-1;



        % lower/upper bounds for parameters
        lb=[ -1 0  0 0 0 0 0 0 0 0 -1 ];
        ub=[ 1  1  1 1 1 1 1 1 1 1 1 ]  ;

        pars_1 = nan(n_iteration,11);
        outp_1 = nan(n_iteration,1);
        for tt=1:n_iteration

            tt

            % random startiing point
            stp=[ leak_range(tt)...
                rand(1,1).*0.9+0.02 ...
                (0.8).*rand(1,1)+0.02...
                (0.9).*rand(1,1)+0.02...
                (0.8).*rand(1,1)+0.02......
                (0.9).*rand(1,1)+0.02 ...
                (0.9).*rand(1,1)+0.02...
                (0.9).*rand(1,1)+0.02 ...
                (0.2).*rand(1,1)+0.15...
                (0.9).*rand(1,1)+0.02...
                (0.9).*rand(1,1)+0.02...
                ];




            % evaluate nll
            [nll]=cost_function_model_9(stp);


            pars_1(tt,:)=stp;
            outp_1(tt)=nll;

        end

        % best starting point parameters
        [k,l]=min(outp_1);
        stp_1=pars_1(l,:);




        plb=stp_1 - abs(0.6*stp_1);
        plb(plb<=lb)=lb(plb<=lb)+0.01;
        pub=stp_1 + abs(0.6*stp_1);
        pub(pub>=ub)=ub(pub>=ub)-0.01;


        options.Display='off';



        options.UncertaintyHandling = true;
        options.NoiseSize = 1;

        % fit
        [xbest,FVAL] = bads('cost_function_model_9',stp_1,[lb ],[ub],[plb],[pub],[],options);








        % test trials
        global dat;
        dat.t=10;


        dat.response=Resp_ses(id_test)';
        dat.right=All_contrasts(id_test,1:10)';
        dat.left=All_contrasts(id_test,11:20)';
        dat.rt=All_RT(id_test);
        dat.conf=All_conf(id_test);
        dat.ses_id=All_ses_id(id_test,:);


        % evaluate the model withh best parameters
        [nll,p,conf,output] = cost_function_model_9(xbest);


        %  best parameters o for this fold

        All_param(cv,:)=xbest;

        %  nll on test trials for this fold
        LL_cv(cv,:)=-nll;

        % concatenate test trials info for this fold
        All_left=[ All_left;dat.left'];
        All_right=[ All_right;dat.right'];
        All_stimulus=[ All_stimulus;dat.right'-dat.left'];
        All_ch_model=[ All_ch_model;p];
        All_Conf=[All_Conf;conf];
        Real_data=[Real_data;Resp_ses(id_test) All_conf(id_test) All_RT(id_test) All_gm(id_test)];


        toc

        % average best parameters across folds


    end
end

% average best parameters across folds

xbest=mean(All_param,1);




% save parameters and fits







modelName = '_model_9';

Parameters = xbest;

save(fullfile(resultsFolder,...
    ['Parameters' modelName '_sbj_' num2str(subject_n) '.mat']),...
    'Parameters')

save(fullfile(resultsFolder,...
    ['All_data' modelName '_sbj_' num2str(subject_n) '.mat']),...
    'All_Conf','All_left','All_right',...
    'All_stimulus','All_ch_model','Real_data')

end


