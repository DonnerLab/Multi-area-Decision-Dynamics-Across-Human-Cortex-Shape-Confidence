function [  Final_DV P_left P_right Choices Confidence output thr_collapse Confidence_end_2 Confidence_left Confidence_right P_cross] = model_simulation(pars,data,bound,n_simulations)

%% simulates N times a model with parameters (pars)
% model types: - 1: no bound, no asymmetruc oocnfidence weight, no FFI
               % -2: fixed bound, FFI, no asymmetry coonf weight
              % -3: collapsing bound, FFI, no asymmetry coonf weight
               % -4: collapsing bound, FFI,  asymmetry coonf weight (best model)
 
%
% n_simulations (optional): number of stochastic simulations per trial
%                           (default: 250)
%%
if nargin<4 || isempty(n_simulations)
    n_simulations=250;
end

tau_a=(pars(3)).*15+0.0001;
Adapt=pars(4)+(1-pars(4))*exp(-(0:9)./tau_a);
%Adapt=data.adapt;
Adapt=Adapt./sum(Adapt);

%
t=data.t-1;

% sigma noise
ext=pars(2)/20;

% add conf asymmetry weight
if bound==4
w_conf=pars(10);
else
w_conf=1;
    
end

% collpasing  bound
dv=abs(pars(6)*50);
thr_1=abs(pars(8)).*15;

if bound==1
    
    thr_collapse=inf(1,10).*thr_1;
    
elseif bound==2
    
    thr_collapse=ones(1,10).*thr_1;
    
else
    
    thr_collapse=(thr_1-(thr_1-(0*thr_1)).*((0:9)./((0:9)+(1./dv))));
    
end

% adapted contrast values

left=data.left(1:t+1,:).*Adapt';
right=data.right(1:t+1,:).*Adapt';


% leak of integration
leak=pars(1);


% iterate thorough trial and get choice for the model
P_left=zeros(10,size(left,2),n_simulations).*left;
P_right=zeros(10,size(right,2),n_simulations).*right;

% mutual inhibition weight
if bound==1
W_FFI=0;

else


W_FFI=pars(5)*1;

end
k=1;



Noise_L=(normrnd(0, ext,[1 size(P_left,2) size(P_left,3)]));
Noise_R=(normrnd(0, ext,[1 size(P_left,2) size(P_left,3)]));

L_to_R=W_FFI.*(left(k,:,:));
R_to_L=W_FFI.*(right(k,:,:));

P_left(k,:,:)=(left(k,:,:)-R_to_L)+Noise_L;
P_right(k,:,:)=(right(k,:,:)-L_to_R)+Noise_R;

% integrate evidence with leak + noise

   

for k=2:t+1
    
    
    Noise_L=(normrnd(0, ext,[1 size(P_left,2) size(P_left,3)]));
    Noise_R=(normrnd(0, ext,[1 size(P_left,2) size(P_left,3)]));
    
   
    
   
 L_to_R=W_FFI.*(left(k,:,:));
R_to_L=W_FFI.*(right(k,:,:));
   
 
    P_left(k,:,:)=(1-leak)*(P_left(k-1,:,:))+(left(k,:,:)-R_to_L)+Noise_L;
    P_right(k,:,:)=(1-leak)*(P_right(k-1,:,:))+(right(k,:,:)-L_to_R)+Noise_R;

   
    
end

% 
B=pars(11);
P_left=P_left+B;
P_right=P_right-B;


P=(P_right-P_left);


P_before=(P(1:end,:,:));

P_end=(P_left(1:end,:,:));
P_abs=cat(4,P_left,P_right);



P_cross=zeros(1,size(P,2),size(P,3));

Ch_end=nan(1,size(P,2),size(P,3));

Ch_L=nan(10,size(P_left,2),size(P_left,3));
id_cross_thr_pos=(P_left>=thr_collapse');
Ch_L(id_cross_thr_pos)=-1;



Ch_R=nan(10,size(P_right,2),size(P_right,3));
id_cross_thr_pos=(P_right>=thr_collapse');
Ch_R(id_cross_thr_pos)=1;


id_both_nan=isnan(Ch_R)&isnan(Ch_L);

Ch=nansum(cat(4,Ch_R,Ch_L),4);

Ch(id_both_nan)=nan;


Confidence_end=nan(1,size(P,2),size(P,3));
Confidence_left=nan(1,size(P,2),size(P,3));
Confidence_right=nan(1,size(P,2),size(P,3));


for tr=1:size(Ch,2)
    
    
    
    nsample=10;
    id_end=find(sum(isnan(Ch(:,tr,:)),1)==10);
    
    
    Ch_end(1,tr,id_end)=sign(P(nsample,tr,id_end));
    P_cross(1,tr,id_end)=1;
    id_z=find(Ch_end(1,tr,id_end)==0);
    Ch_end(1,tr,id_end(id_z))=sign(randn(size(id_z)));
    
    id_pos=find(Ch_end(1,tr,id_end)==1);
    Confidence_end(1,tr,id_end(id_pos))=abs(P_right(nsample,tr,id_end(id_pos))-w_conf.*P_left(nsample,tr,id_end(id_pos)));
    id_neg=find(Ch_end(1,tr,id_end)==-1);
    Confidence_end(1,tr,id_end(id_neg))=abs(P_left(nsample,tr,id_end(id_neg))-w_conf.*P_right(nsample,tr,id_end(id_neg)));
    
    
    id_pos=find(Ch_end(1,tr,id_end)==1);
    Confidence_left(1,tr,id_end(id_pos))=(P_left(nsample,tr,id_end(id_pos)));
    id_neg=find(Ch_end(1,tr,id_end)==-1);
    Confidence_left(1,tr,id_end(id_neg))=(P_left(nsample,tr,id_end(id_neg)));
    
    id_pos=find(Ch_end(1,tr,id_end)==1);
    Confidence_right(1,tr,id_end(id_pos))=(P_right(nsample,tr,id_end(id_pos)));
    id_neg=find(Ch_end(1,tr,id_end)==-1);
    Confidence_right(1,tr,id_end(id_neg))=(P_right(nsample,tr,id_end(id_neg)));
    
    
end


P_cross=zeros(size(Ch));

%%
for i=1:10
    
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
    
    id_pos=find(this_Choice==1);
    Confidence_left(id_pos)=(this_P_left(id_pos));
    id_neg=find(this_Choice==-1);
    Confidence_left(id_neg)=(thr_collapse(i));
    
    
    id_pos=find(this_Choice==1);
    Confidence_right(id_pos)=(thr_collapse(i));
    id_neg=find(this_Choice==-1);
    Confidence_right(id_neg)=(this_P_right(id_neg));

    
   prev_P_right=P_right(i,:,:);
    prev_P_left=P_left(i,:,:);   
    
      id_pos_eq=find(this_Choice(id_equal)==1);
     Confidence_end(id_equal(id_pos_eq))=((this_P_right(id_equal(id_pos_eq)))-w_conf.*(this_P_left(id_equal(id_pos_eq))));     
     Confidence_right(id_equal)=(this_P_right(id_equal));
     Confidence_left(id_equal)=(this_P_left(id_equal));

     
     id_neg_eq=find(this_Choice(id_equal)==-1);
     Confidence_end(id_equal(id_neg_eq))=((this_P_left(id_equal(id_neg_eq)))-w_conf.*(this_P_right(id_equal(id_neg_eq))));     
 
    
    id=(Ch(i,:,:)==1|Ch(i,:,:)==-1|Ch(i,:,:)==0);
    Ch(i+1:end,id)=nan;
    P_cross(i,id)=1;
    P_cross(i+1:end,id)=0;

end


%%
Ch_end(Ch_end==-1)=0;


slope_conf=abs(pars(7)*30);
thr_conf=abs(pars(9)*4);   


Conf_all=(abs(Confidence_end(:,:,:)));
Final_DV=Conf_all;

[min_v]=min(Conf_all,[],2);min_v=min(min_v,[],1);
[max_v]=max(Conf_all,[],2);max_v=max(max_v,[],1);

min_v=repmat(min_v,1,size(Conf_all,2),1);
max_v=repmat(max_v,1,size(Conf_all,2),1);

Norm_Conf_all=((Conf_all-min_v)./(max_v-min_v));



Median=median(Norm_Conf_all,2);
Confidence_end_2=Norm_Conf_all;

Confidence_end(:,:,:)=1./(1+exp(-slope_conf.*((Norm_Conf_all-Median.*thr_conf))));
id_equal=find((abs(Confidence_end(:,:,:))==0.5));
Confidence_end=double(Confidence_end>0.5);
Confidence_end(id_equal)=randi(2,size(id_equal))-1;





id_1_1=(Ch_end(1,:,:)==1&Confidence_end(1,:,:)==1);
P_confidence(:,1)=sum(id_1_1,3)./n_simulations;
id_1_0=(Ch_end(1,:,:)==1&Confidence_end(1,:,:)==0);
P_confidence(:,2)=sum(id_1_0,3)./n_simulations;
id_0_1=(Ch_end(1,:,:)==0&Confidence_end(1,:,:)==1);
P_confidence(:,3)=sum(id_0_1,3)./n_simulations;
id_0_0=(Ch_end(1,:,:)==0&Confidence_end(1,:,:)==0);
P_confidence(:,4)=sum(id_0_0,3)./n_simulations;



Choices=Ch_end;
output=Confidence_end;

Confidence=P_confidence;
%%
end


