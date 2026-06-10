function Decoder_Performance = decode_choice_all_betas(subject,sessions,Rec)


Extra_sessions=[1,2,4,5,7,0,11,12,15,17,18,10,0,0,0,0,0,0,0,0,0];

% if session 1 was re-an use extra session
if sessions==1 & sum(Extra_sessions==subject)==1
    sessions=11;
end

% areas from which to decode
area={'4','6a','6d'};

agg_dir='~/Desktop/2AFC/AGG/400/matlab_pairs/';
behav_dir='~/Desktop/2AFC/trials/';

addpath('~/Desktop/2AFC/Analyses/liblinear-2.45/matlab')
addpath('~/aletoso/Desktop/2AFC/trials/')

if subject<10
    subject=['S0',num2str(subject)];
else
    subject=['S',num2str(subject)];
end


TFR_area_L_all=[];
All_ERF_R=[];
All_Behavior=[];
Samples=[];


sides={'l','r'};

%  load relevant ERF and Behavioral data

for ses=1:size(sessions,2)
    for r=1:size(Rec,2)
        
        cd('~/Desktop/2AFC/ERF/')
        
        path=['~/Desktop/2AFC/ERF/',subject,'SES',num2str(sessions(ses)),'_',num2str(Rec(r)),'_',cell2mat(area(1)),'_ERF_r.mat'];
        
        if isfile(path)
            
            All_ERF=[];
            for side=1:2
                
                %load ERF from subject s and concatenate across areas and
                %sessions
                
                
                cd('~/Desktop/2AFC/ERF/')
                
                for i=1:size(area,2)
                    path=['~/Desktop/2AFC/ERF/',subject,'SES',num2str(sessions(ses)),'_',num2str(Rec(r)),'_',cell2mat(area(i)),'_ERF_',cell2mat(sides(side)),'.mat'];
                    if isfile(path)
                        
                        load(['~/Desktop/2AFC/ERF/',subject,'SES',num2str(sessions(ses)),'_',num2str(Rec(r)),'_',cell2mat(area(i)),'_ERF_',cell2mat(sides(side)),'.mat']);
                        % erf is trialsXvertexXtime
                        erf=(comb_dict{1,1}.erfdata(:,:,:));
                        
                        
                        if isempty(erf)
                            
                        else
                            All_ERF=cat(2,erf,All_ERF);
                            
                        end
                        
                        
                    else
                    end
                    
                    
                end
                % save N of vertices per hemifield
                n_vx_side(side,:)=size(All_ERF,2);
            end
            
            % load behavioral results
            cd('~/Desktop/2AFC/trials/')
            if isfile(path)
                    load(['~/Desktop/2AFC/trials/',subject,'_',num2str(sessions(ses)),'_',num2str(Rec(r)),'_trials.mat']);
                
            else
            end
            
            data.trial=All_ERF;
            data.label=cellstr(num2str([1:size(All_ERF,2)]'));
            
            All_ERF=permute(All_ERF,[3,1,2]);
            All_ERF=permute(All_ERF,[2,3,1]);
            
            % baseline from pre-stim interval
            All_ERF_R=cat(1,All_ERF_R,abs((All_ERF)-(mean(All_ERF(:,:,1:201),3))));
            All_Behavior=[All_Behavior;data.trialinfo];
            
            
            
            
            
        else
        end
    end
end



tic



    
    
        % extract trial index for 10-fold cv
        id_cv=cvpartition(size(All_Behavior,1),'Kfold',10);
        
        % train and test decoder separately for each time t
        for t=1:size(All_ERF_R,3)
            
            X_both=[];Y=[];
            
            % source-reconstructed ERF
            X=[All_ERF_R(:,:,t)];
            
            
            X_all=X;
            
                % choice           
            Y_all=double(All_Behavior(:,2));
            Y_all(Y_all==0)=-1;
            
            yFit_all=[];real_Y=[];Perf_all=[];
            Beta_all=[];All_samples=[];All_Choices=[];block_type=[];
            
            % 10 folded cv
            for cv=1:10
                
                Y=[];X=[];
                
                % extract training/test trials for each cv fold
                id_train=training(id_cv,cv);
                id_test=test(id_cv,cv);
                
                Y=Y_all(id_train);
                X=X_all(id_train,:);
                
                % up-sample minority choice category in traning data
                id_left=find(Y==-1);id_right=find(Y==1);
                if size(id_left,1)>size(id_right,1)
                    n_extra=size(id_left,1)-size(id_right,1);
                    id_extra=id_right(randi(size(id_right,1),n_extra,1));
                else
                    n_extra=size(id_right,1)-size(id_left,1);
                    id_extra=id_left(randi(size(id_left,1),n_extra,1));
                end
                
                Y=[Y;Y(id_extra)];
                X=[X;X(id_extra,:)];
                
                
                % z-score training X data
                [Z,mu,sigma]=zscore(X);
                
                % train logistic model
                Mdl=fitclinear(zscore(X),(Y),'Learner','logistic');%'logiRegularization','ridge','Lambda',[1]);
                
                
                % predict test choices (z-score test X data using mu/sigma from train data)
                Y_fit=predict(Mdl,normalize(X_all(id_test,:),'center',mu,'scale',sigma));
                
                % extract each vx beta           
                Beta_fit=Mdl.Beta;

                % concatenate alcross cv-fold real test Y, predicted Y,
                % prob of choice==1, Stimulus ID, block type
                yFit_cv=Y_fit;
                test_y=[];
                test_y=(Y_all(id_test));
                yFit_all=[yFit_all;yFit_cv];
                
                real_Y=[real_Y;test_y];
                
                Beta_all=[Beta_all;Beta_fit'];
                
                All_samples=[All_samples;All_Behavior(id_test,16:35)];
                All_Choices=[All_Choices;All_Behavior(id_test,2)];
                block_type=[block_type;All_Behavior(id_test,13)];
                
            end
            
            
            % Calculate performance
            Perfs=mean((yFit_all)==(real_Y));
            
            Beta(1,t,:)=mean(Beta_all,1);
            
            Decoder_Performance_bt(1,t,1)=Perfs(1);
        end
    



% save

Decoder_Performance=(Decoder_Performance_bt);

cd   ~/Desktop/2AFC/Analyses/Decoder_Final/ERF/PCA/all_Vx
name=[(subject),'_train_CHOICE_decoder_session',num2str(sessions),'_PCA_PMd_M1_SINGLETRIAL.mat'];

save(name,'Beta','All_samples','All_Choices','block_type','All_Behavior','n_vx_side')

cd   ~/Desktop/2AFC/Analyses/Decoder_Final/ERF/PCA/all_Vx
name=[(subject),'_CHOICE_decoder_session',num2str(sessions),'_PCA_PMd_M1.mat'];

save(name,'Decoder_Performance')

toc
delete(gcp('nocreate'));

end