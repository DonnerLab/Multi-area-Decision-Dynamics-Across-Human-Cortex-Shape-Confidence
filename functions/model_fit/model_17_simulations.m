function [ Choices Confidence output] = model_17_simulations(pars,data)
%%
% number of simulations
n_simulations=1000;

% adaptation kernel        
tau_a=(pars(3)).*15+0.0001;
Adapt=pars(4)+(1-pars(4))*exp(-(0:9)./tau_a);
Adapt=Adapt./sum(Adapt);

% 
t=data.t-1;

% sigma noise
ext=pars(2)/20;


% collapsing  bound
dv=abs(pars(6)*50);
thr_1=abs(pars(8)).*15;
thr_collapse=(thr_1-(thr_1-(0*thr_1)).*((0:9)./((0:9)+(1./dv))));


% apply adaptation kernels
 left=data.left(1:t+1,:).*Adapt';
 right=data.right(1:t+1,:).*Adapt';
 
% leak of integration
leak=pars(1);


% iterate thorough trial and get choice for the model
P_left=zeros(10,size(left,2),n_simulations).*left;
P_right=zeros(10,size(right,2),n_simulations).*right;

% inhibition  weight

w_I=pars(5)*1;

k=1;
 

 
  Noise_L=(normrnd(0, ext,[1 size(P_left,2) size(P_left,3)]));
    Noise_R=(normrnd(0, ext,[1 size(P_left,2) size(P_left,3)]));
    
    % side parameer
B=pars(11).*1;
   
    
% first sample is accumulated

P_left(k,:,:)=(left(k,:,:))+Noise_L;
P_right(k,:,:)=(right(k,:,:))+Noise_R;

 % 0-rectification

     P_left(P_left<0)=0;
     P_right(P_right<0)=0;

% integrate evidence with leak + noise
for k=2:t+1
    
 % left and right noise
  Noise_L=(normrnd(0, ext,[1 size(P_left,2) size(P_left,3)]));
  Noise_R=(normrnd(0, ext,[1 size(P_left,2) size(P_left,3)]));
    
  % lca
L_to_R=w_I.*(P_left(k-1,:,:));
R_to_L=w_I.*(P_right(k-1,:,:));


   
    % integration
    P_left(k,:,:)=(1-leak)*(P_left(k-1,:,:))+(left(k,:,:)-R_to_L)+Noise_L;
    P_right(k,:,:)=(1-leak)*(P_right(k-1,:,:))+(right(k,:,:)-L_to_R)+Noise_R;


     % 0-rectification

     P_left(P_left<0)=0;
     P_right(P_right<0)=0;

end  

% add fixed bias
P_left=P_left+B;
P_right=P_right-B;



P=(P_right-P_left);




      
 Ch_end=nan(1,size(P,2),size(P,3));   
  
 % find bound crossed left DV
 Ch_L=nan(10,size(P_left,2),size(P_left,3));   
 id_cross_thr_pos=(P_left>=thr_collapse');
 Ch_L(id_cross_thr_pos)=-1;

 % find bound crossed right DV
 Ch_R=nan(10,size(P_right,2),size(P_right,3));   
 id_cross_thr_pos=(P_right>=thr_collapse');
 Ch_R(id_cross_thr_pos)=1;

% no crossing in both DV
 id_both_nan=isnan(Ch_R)&isnan(Ch_L); 
 Ch=nansum(cat(4,Ch_R,Ch_L),4);   
 Ch(id_both_nan)=nan;

Confidence_end=nan(1,size(P,2),size(P,3));    

% trials with no bound cross
for tr=1:size(Ch,2)

   
   nsample=10;
    
 id_end=find(sum(isnan(Ch(:,tr,:)),1)==10);

 % choice is the higher DV
Ch_end(1,tr,id_end)=sign(P(nsample,tr,id_end));

% if zero difference set choice randomly
id_z=find(Ch_end(1,tr,id_end)==0);
Ch_end(1,tr,id_end(id_z))=sign(randn(size(id_z)));

% confidence value at stimulus offset
id_pos=find(Ch_end(1,tr,id_end)==1);
Confidence_end(1,tr,id_end(id_pos))=abs(P_right(nsample,tr,id_end(id_pos))-P_left(nsample,tr,id_end(id_pos)));
id_neg=find(Ch_end(1,tr,id_end)==-1);
Confidence_end(1,tr,id_end(id_neg))=abs(P_left(nsample,tr,id_end(id_neg))-P_right(nsample,tr,id_end(id_neg)));
 
end
  
% trials with  bound crossing
for i=1:10
   %%
% check at each sample n when bound is crossed
    
    this_Choice=Ch(i,:,:);
    this_P_right=P_right(i,:,:);
    this_P_left=P_left(i,:,:);

    id=find(this_Choice==1|this_Choice==-1);
    Ch_end(id)=this_Choice(id);
    
    id_equal=find(this_Choice==0);    
    Ch_end(id_equal)=double(this_P_right(id_equal)>this_P_left(id_equal));  
    
    this_Choice(id_equal)=Ch_end(id_equal);
    id=find(this_Choice==0);
    this_Choice(id)=-1;
    
    
   
  
    
    
    Ch_end(Ch_end==0)=-1;
   
    id_pos=find(this_Choice==1);
    Confidence_end(id_pos)=thr_collapse(i)-w_conf.*(this_P_left(id_pos));
    id_neg=find(this_Choice==-1);
    Confidence_end(id_neg)=(thr_collapse(i)-w_conf.*(this_P_right(id_neg)));
    
     id_pos_eq=find(this_Choice(id_equal)==1);
     Confidence_end(id_equal(id_pos_eq))=((this_P_right(id_equal(id_pos_eq)))-(this_P_left(id_equal(id_pos_eq))));     
    id_neg_eq=find(this_Choice(id_equal)==-1);
     Confidence_end(id_equal(id_neg_eq))=((this_P_left(id_equal(id_neg_eq)))-(this_P_right(id_equal(id_neg_eq))));     
 
    
    id=(Ch(i,:,:)==1|Ch(i,:,:)==-1|Ch(i,:,:)==0);
    Ch(i+1:end,id)=nan;
    
end

Ch_end(Ch_end==-1)=0;

% transform confidence in binary confidence reports
slope_conf=abs(pars(7)*30);
thr_conf=abs(pars(9)*4);
    
% normalize confidence from 0 to 1
 Conf_all=(abs(Confidence_end(:,:,:)));
 [min_v]=min(Conf_all,[],2);min_v=min(min_v,[],1);
 [max_v]=max(Conf_all,[],2);max_v=max(max_v,[],1);

 min_v=repmat(min_v,1,size(Conf_all,2),1);
 max_v=repmat(max_v,1,size(Conf_all,2),1);

 Norm_Conf_all=((Conf_all-min_v)./(max_v-min_v));
   

 % logistic function
 Median=median(Norm_Conf_all,2);
      
Confidence_end(:,:,:)=1./(1+exp(-slope_conf.*((Norm_Conf_all-Median.*thr_conf))));

% binarize
id_equal=find((abs(Confidence_end(:,:,:))==0.5));
Confidence_end=double(Confidence_end>0.5);
Confidence_end(id_equal)=randi(2,size(id_equal))-1;


% probability of each confidence and choice category
   id_1_1=(Ch_end(1,:,:)==1&Confidence_end(1,:,:)==1);
    P_confidence(:,1)=(sum(id_1_1,3))./(n_simulations);
    id_1_0=(Ch_end(1,:,:)==1&Confidence_end(1,:,:)==0);
    P_confidence(:,2)=(sum(id_1_0,3)./n_simulations);
    id_0_1=(Ch_end(1,:,:)==0&Confidence_end(1,:,:)==1);
    P_confidence(:,3)=(sum(id_0_1,3)./n_simulations);
    id_0_0=(Ch_end(1,:,:)==0&Confidence_end(1,:,:)==0);
    P_confidence(:,4)=(sum(id_0_0,3)./n_simulations);


% P choice R>L

Choices=mean(Ch_end,3)';
output=sum([P_confidence(:,1),P_confidence(:,3)],2)';

Confidence=P_confidence;

end