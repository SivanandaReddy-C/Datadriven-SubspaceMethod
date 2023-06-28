%% This code is on Subspace and Difference method frameworks for Motor Imagery EEG classification.
% Paper: Data-driven Motor Imagery EEG Classifier using Difference Subspace Method
% Authors: C Sivananda Reddy, M Ramasubba Reddy.
% Please contact us if you find any issues with this code.
% Email: sivananda.reddi@gmail.com


% i. It is implemented on BCI competition 4 dataset 2b. 
% ii. It is a generalized code which can be implemented on any datasets. 
% iii.The dimesion of the subspaces is tried on all values, which can be 
%     replaced with the optimum one by working on suitable selection algorithms

clear all;
close all;
clc;

% BCI 4_2B dataset
fs=250;
addpath([pwd '/../dataset/BCI4_2B/']);

[s_sess1,h_sess1]=sload('B0101T.gdf'); % B0101T dataset preparation
Trials_All_Ind_sess1=find(h_sess1.EVENT.TYP==768);
Reject_Ind_sess1=find(h_sess1.ArtifactSelection==1);

Trials_All_Ind_sess1(Reject_Ind_sess1)=[];
h_sess1.Classlabel(Reject_Ind_sess1)=[];

Trials_Ind_sess1=Trials_All_Ind_sess1;
Z_sess1= h_sess1.Classlabel';
trial=1;
for i=1:length(Trials_Ind_sess1)
    st_sess1=h_sess1.EVENT.POS(Trials_Ind_sess1(i));
    dur_sess1=h_sess1.EVENT.DUR(Trials_Ind_sess1(i));
    Y_temp_sess1(:,:)=s_sess1(st_sess1:st_sess1+dur_sess1-1,1:3);
    if(find(isnan(Y_temp_sess1(:,:)))~=0) 
       i
       Z_sess1(trial)=[];
    else
        Y_sess1_mu(:,:,trial)=(bandpass(Y_temp_sess1(:,:),[8 12],fs))';
        Y_sess1_beta(:,:,trial)=(bandpass(Y_temp_sess1(:,:),[12 30],fs))';
        trial=trial+1;
    end
end

[s_sess2,h_sess2]=sload('B0102T.gdf'); % B0102T dataset preparation
Trials_All_Ind_sess2=find(h_sess2.EVENT.TYP==768);
Reject_Ind_sess2=find(h_sess2.ArtifactSelection==1);

Trials_All_Ind_sess2(Reject_Ind_sess2)=[];
h_sess2.Classlabel(Reject_Ind_sess2)=[];

Trials_Ind_sess2=Trials_All_Ind_sess2;
Z_sess2= h_sess2.Classlabel';
trial=1;
for i=1:length(Trials_Ind_sess2)
    st_sess2=h_sess2.EVENT.POS(Trials_Ind_sess2(i));
    dur_sess2=h_sess2.EVENT.DUR(Trials_Ind_sess2(i));
    Y_temp_sess2(:,:)=s_sess2(st_sess2:st_sess2+dur_sess2-1,1:3);
    if(find(isnan(Y_temp_sess2(:,:)))~=0) 
       i
       Z_sess2(trial)=[];
    else
        Y_sess2_mu(:,:,trial)=(bandpass(Y_temp_sess2(:,:),[8 12],fs))';
        Y_sess2_beta(:,:,trial)=(bandpass(Y_temp_sess2(:,:),[12 30],fs))';
        trial=trial+1;
    end
end

[s_sess3,h_sess3]=sload('B0103T.gdf'); % B0103T dataset preparation
Trials_All_Ind_sess3=find(h_sess3.EVENT.TYP==768);
Reject_Ind_sess3=find(h_sess3.ArtifactSelection==1);

Trials_All_Ind_sess3(Reject_Ind_sess3)=[];
h_sess3.Classlabel(Reject_Ind_sess3)=[];

Trials_Ind_sess3=Trials_All_Ind_sess3;
Z_sess3= h_sess3.Classlabel';
trial=1;
for i=1:length(Trials_Ind_sess3)
    st_sess3=h_sess3.EVENT.POS(Trials_Ind_sess3(i));
    dur_sess3=h_sess3.EVENT.DUR(Trials_Ind_sess3(i));
    Y_temp_sess3(:,:)=s_sess3(st_sess3:st_sess3+1999,1:3);
    if(find(isnan(Y_temp_sess3(:,:)))~=0) 
       i
       Z_sess3(trial)=[];
    else
        Y_sess3_mu(:,:,trial)=(bandpass(Y_temp_sess3(:,:),[8 12],fs))';
        Y_sess3_beta(:,:,trial)=(bandpass(Y_temp_sess3(:,:),[12 30],fs))';
        trial=trial+1;
    end
end

[s_sess4,h_sess4]=sload('B0104E.gdf'); % B0104T dataset preparation
B0404E_Labels=load('B0104E.mat');

Trials_All_Ind_sess4=find(h_sess4.EVENT.TYP==768);
Reject_Ind_sess4=find(h_sess4.ArtifactSelection==1);

Trials_All_Ind_sess4(Reject_Ind_sess4)=[];
B0404E_Labels.classlabel(Reject_Ind_sess4)=[];

Trials_Ind_sess4=Trials_All_Ind_sess4;
Z_sess4= B0404E_Labels.classlabel';
trial=1;
for i=1:length(Trials_Ind_sess4)
    clear Y_temp;
    st_sess4=h_sess4.EVENT.POS(Trials_Ind_sess4(i));
    dur_sess4=h_sess4.EVENT.DUR(Trials_Ind_sess4(i));
    Y_temp_sess4(:,:)=s_sess4(st_sess4:st_sess4+1999,1:3);
    if(find(isnan(Y_temp_sess4(:,:)))~=0) 
       i
       Z_sess4(trial)=[];
    else
        Y_sess4_mu(:,:,trial)=(bandpass(Y_temp_sess4(:,:),[8 12],fs))';
        Y_sess4_beta(:,:,trial)=(bandpass(Y_temp_sess4(:,:),[12 30],fs))';
        trial=trial+1;
    end
end

[s_sess5,h_sess5]=sload('B0105E.gdf'); % B0403T dataset preparation
B0405E_Labels=load('B0105E.mat');

Trials_All_Ind_sess5=find(h_sess5.EVENT.TYP==768);
Reject_Ind_sess5=find(h_sess5.ArtifactSelection==1);

Trials_All_Ind_sess5(Reject_Ind_sess5)=[];
B0405E_Labels.classlabel(Reject_Ind_sess5)=[];

Trials_Ind_sess5=Trials_All_Ind_sess5;
Z_sess5= B0405E_Labels.classlabel';
trial=1;
for i=1:length(Trials_Ind_sess5) 
    clear Y_temp;
    st_sess5=h_sess5.EVENT.POS(Trials_Ind_sess5(i));
    dur_sess5=h_sess5.EVENT.DUR(Trials_Ind_sess5(i));
    Y_temp_sess5(:,:)=s_sess5(st_sess5:st_sess5+1999,1:3);
    if(find(isnan(Y_temp_sess5(:,:)))~=0) 
       i
       Z_sess5(trial)=[];
    else
        Y_sess5_mu(:,:,trial)=(bandpass(Y_temp_sess5(:,:),[8 12],fs))';
        Y_sess5_beta(:,:,trial)=(bandpass(Y_temp_sess5(:,:),[12 30],fs))';
        trial=trial+1;
    end
end

st=3.5*fs+1:1:5*fs+1; 
ed=5.5*fs:1:7*fs; % Start and end samples of 2sec window
Ch=3;

% Training data: First 3 sessions, Test data: Session 4.
YTn=cat(3,[Y_sess1_mu;Y_sess1_beta],[Y_sess2_mu;Y_sess2_beta],[Y_sess3_mu;Y_sess3_beta]);
YTt=[Y_sess4_mu;Y_sess4_beta];
ZTn=[Z_sess1 Z_sess2 Z_sess3];
ZTt=Z_sess4;

indL=find(ZTn==1); % Indices of class1 and class2 data
indR=find(ZTn==2);
NL=length(indL);
NR=length(indR);
L=min([NL,NR]);
indTt_L=find(ZTt==1);
indTt_R=find(ZTt==2);
NTest=length(ZTt);

%% i.Subspace method 
NRow=1;
for d=1:2*Ch-1
    clear W;
    for epoch=1:length(st)
        st_time=st(epoch);
        ed_time=ed(epoch);
    
        % Covariance computation of class 1 data
        clear YL; clear Rl; 
        Rl=0;
        for trial=1:L 
            Ytemp(:,:)=YTn(:,st_time:ed_time,indL(trial));
            Rl=Rl+cov(Ytemp');
        end
        Rl=Rl./L;
    
        % Covariance computation of class 2 data
        clear YR; clear Rr;
        Rr=0;
        for trial=1:L 
            Ytemp(:,:)=YTn(:,st_time:ed_time,indR(trial));
            Rr=Rr+cov(Ytemp');
        end
        Rr=Rr./L;
    
        % Bases computation of class 1 and class 2 subspaces
        [Ur Er ~]=svd(Rr); [Ul El ~]=svd(Rl); 
    
        % Measuring orthogonality between class1 and class2 subspaces
        clear Url; clear Erl; clear Srl;
        [Url Erl ~]=svd(Ur(:,1:d)'*Ul(:,1:d)); cos_theta=diag(Erl);
        Srl=sum(cos_theta(1:d).^2)/d;
        W(epoch)=1-(Srl); 
    end
    
    [VV II]=max(W);
    st_time=st(II);
    ed_time=ed(II); % 2seconds segment of maximum orthogonality between class1 and class2 subspaces
    
    % Linear Subspace model of optimum 2seconds segment
    % Covariance computation 
    clear YL; clear Rl; 
    Rl=0;
    for trial=1:L 
        Ytemp(:,:)=YTn(:,st_time:ed_time,indL(trial));
        Rl=Rl+cov(Ytemp');
    end
    Rl=Rl./L;
    
    clear YR; clear Rr;
    Rr=0;
    for trial=1:L 
        Ytemp(:,:)=YTn(:,st_time:ed_time,indR(trial));
        Rr=Rr+cov(Ytemp');
    end
    Rr=Rr./L;
    
    % Optimum bases for the class1 and class2 subspaces 
    [Ur Er ~]=svd(Rr); [Ul El ~]=svd(Rl); 
    
    % Classifying test trial based on its closeness with the class1 and class2
    % subspaces
    for epoch=1:length(st)
       st_time=st(epoch);
       ed_time=ed(epoch);
       clear Zpredict;
       for trial=1:NTest
            Ytemp(:,:)=YTt(:,st_time:ed_time,trial);
            RTt=Ytemp*Ytemp';
            [UTt ETt ~]=svd(RTt);
            clear cos_theta; clear ETt_l; clear STt_l;
            [~,ETt_l,~]=svd(UTt(:,1:d)'*Ul(:,1:d));    cos_theta=diag(ETt_l);
            STt_l=sum(cos_theta(1:d).^2)/d;
            clear cos_theta; clear ETt_r; clear STt_r;
            [~,ETt_r,~]=svd(UTt(:,1:d)'*Ur(:,1:d));    cos_theta=diag(ETt_r);
            STt_r=sum(cos_theta(1:d).^2)/d;
            [VV,Zpredict(trial)]=max([STt_l,STt_r]);
       end
       [Cf order]=confusionmat(ZTt,Zpredict');
       Acc(epoch)=(Cf(1,1)+Cf(2,2))/NTest;  
       Cf_mat(:,:,epoch)=Cf;
    end
    [Vend Iend]=max(Acc);
    stats=statsOfMeasure(Cf_mat(:,:,Iend),1)
    Result(NRow,:)=[d max(Acc)]
    NRow=NRow+1;
end

%% ii. Difference subspace method
clear Result;
NROW=1;
for d=1:2*Ch-1
clear WDS;
for epoch=1:length(st)
    st_time=st(epoch);
    ed_time=ed(epoch);

    % Covariance computation of class 1
    clear Rl; clear Ytemp;
    Rl=0;
    for trial=1:L 
        Ytemp(:,:)=YTn(:,st_time:ed_time,indL(trial));
        Rl=Rl+cov(Ytemp');
    end
    Rl=Rl./L;

    % Covariance computation of class 2
    clear Rr; clear Ytemp;
    Rr=0;
    for trial=1:L 
        Ytemp(:,:)=YTn(:,st_time:ed_time,indR(trial));
        Rr=Rr+cov(Ytemp');
    end
    Rr=Rr./L;

    % Bases computation of class1 and class2 subspaces 
    clear Ur; clear Ul; clear Er; clear El;
    [Ur Er ~]=svd(Rr); [Ul El ~]=svd(Rl); 

    % Projections matrices of class 1 and class 2 
    P=0;Q=0;
    for i=1:d
       P=P+Ur(:,i)*Ur(:,i)';
       Q=Q+Ul(:,i)*Ul(:,i)';
    end

    % Bases for sum of subspaces
    clear UDS; clear EDS;
    [UDS,EDS,~]=svd(P+Q);
   
    D=UDS(:,find(diag(EDS)<1));
    clear Ytemp; clear YLDS;  clear RlDS;    
    RlDS=0;
    for trial=1:L  % projecting left data onto DS
        Ytemp(:,:)=YTn(:,st_time:ed_time,indL(trial));
        YLDS=D'*Ytemp;
        RlDS=RlDS+cov(YLDS');
    end
    RlDS=RlDS./L; 

    clear Ytemp; clear YRDS; clear RrDS; 
    RrDS=0;
    for trial=1:L % projecting right data onto DS
         Ytemp(:,:)=YTn(:,st_time:ed_time,indR(trial));
         YRDS=D'*Ytemp;
         RrDS=RrDS+cov(YRDS');
    end
    RrDS=RrDS./L;

    clear UrDS; clear ErDS; clear UlDS; clear ElDS;
    % Right and Left subspaces on diff. subspace
    [UrDS ErDS ~]=svd(RrDS);   [UlDS ElDS ~]=svd(RlDS); 

    % Measuring orthogonality among classes of training data     
     for dDs=1:size(D,2)-1
          clear cos_theta; clear ErlDS; 
          [~,ErlDS,~]=svd(UrDS(:,1:dDs)'*UlDS(:,1:dDs));  cos_theta=diag(ErlDS);
          SrlDS=sum(cos_theta(1:dDs).^2)/dDs;
          WDS(dDs,epoch)=1-SrlDS;
     end
end

% Optimum model and its testing
for dDs=1:size(D,2)-1
     [VV II]=max(WDS(dDs,:)); % Getting optimum epoch
     st_time=st(II);
     ed_time=ed(II);
     % Optimum model
     clear Rl; clear Ytemp;
     Rl=0;
     for trial=1:L 
          Ytemp(:,:)=YTn(:,st_time:ed_time,indL(trial));
          Rl=Rl+cov(Ytemp');
     end
     Rl=Rl./L;

     clear Rr; clear Ytemp;
     Rr=0;
     for trial=1:L 
          Ytemp(:,:)=YTn(:,st_time:ed_time,indR(trial));
          Rr=Rr+cov(Ytemp');
     end
     Rr=Rr./L;

     clear Ur; clear Ul; clear Er; clear El;
     [Ur Er ~]=svd(Rr); [Ul El ~]=svd(Rl); 

     P=0;Q=0; % Difference subspace finding on optimum epoch
     for i=1:d
         P=P+Ur(:,i)*Ur(:,i)';
         Q=Q+Ul(:,i)*Ul(:,i)';
     end
     clear UDS; clear EDS;
     [UDS,EDS,~]=svd(P+Q);

     clear D;     
     D=UDS(:,find(diag(EDS)<1));

     clear Ytemp; clear YLDS; clear RlDS;
     RlDS=0;
     for trial=1:L % projecting left data onto DS
          Ytemp(:,:)=YTn(:,st_time:ed_time,indL(trial));
          YLDS=D'*Ytemp;
          RlDS=RlDS+cov(YLDS');
     end
     RlDS=RlDS./L;

     clear Ytemp; clear YRDS;clear RrDS;
     RrDS=0;
     for trial=1:L % projecting right data onto DS
          Ytemp(:,:)=YTn(:,st_time:ed_time,indR(trial));
          YRDS=D'*Ytemp;
          RrDS=RrDS+cov(YRDS');
     end
     RrDS=RrDS./L;

     % Optimum right and left subspaces on diff. subspace
     clear UrDS; clear UlDS; clear ElDS; clear ErDS;
     [UrDS ErDS ~]=svd(RrDS);   [UlDS ElDS ~]=svd(RlDS);

     % Testing
    for epoch=1:length(st)
         st_time=st(epoch);
         ed_time=ed(epoch);
         for trial=1:NTest
              clear Ytemp; clear YTestDS;
              Ytemp(:,:)=YTt(:,st_time:ed_time,trial);
              YTestDS=D'*Ytemp; % projecting test trial onto DS

              clear RTtDS; clear UTtDS; clear ETtDS;
              RTtDS=YTestDS*YTestDS';
              [UTtDS ETtDS ~]=svd(RTtDS); % Finding subspace of projected test trial

              % Similarity check between test subspace and left subspace on DS
              clear cos_theta; clear ETt_lDS; 
              [~,ETt_lDS,~]=svd(UTtDS(:,1:dDs)'*UlDS(:,1:dDs));    cos_theta=diag(ETt_lDS);
              STt_lDS=sum(cos_theta(1:dDs).^2)/dDs;

              % Similarity check between test subspace and right subspace on DS
              clear cos_theta; clear ETt_rDS;
              [~,ETt_rDS,~]=svd(UTtDS(:,1:dDs)'*UrDS(:,1:dDs));    cos_theta=diag(ETt_rDS);
              STt_rDS=sum(cos_theta(1:dDs).^2)/dDs;
              [VV,Zpredict(trial)]=max([STt_lDS,STt_rDS]);
          end
          [Cf order]=confusionmat(ZTt,Zpredict);
          Acc(epoch)=(Cf(1,1)+Cf(2,2))/NTest;
          Cf_mat(:,:,epoch)=Cf;
    end
    [Vend Iend]=max(Acc);
    stats=statsOfMeasure(Cf_mat(:,:,Iend),1)
    Result(NROW,:)=[d dDs max(Acc)]
    NROW=NROW+1; 
    end
end       


