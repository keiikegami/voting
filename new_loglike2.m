

function [loglik]=new_loglike(param,data1)



%tic
global X;       %data of demographic characteristcs
global DATA;
global dFX      % (joint distribution of dFXRace,dFXEduc,dFXIncm)
global Votes;   % (vote share: clark, dean, edwards, kerry)
global VOther;  %Votes of Penna & Sharpton etc. % こいつら候補者か
global RDemHat; % (registred number of democrats predicted)
global PopTot;  % (total population)
global Cand;    % (candidates on ballot 35 x 4)
global Open;    % (dummy for open primary)
global MOpen;   % (dummy for modified open primary)
global time     %
global n_eval
global ii iii jj jjj kk kkk
global N_dFX    % Number of categories for voter type
global N_sim
global bandwidth
global SimAlp
global SimXsi
global i_bayes
global bestpara
global T

if i_bayes==1
    param=param;
end

% composite (75x1)  : For each state, we have at most 4 values (total 75)
% Tij (111x1)        : For each state, we have at most 6 values (total 111)
% Cx    ((sizeXk-1) x1)  : Cost function parameter.  Number of Xk
% Cz    (3 x 1)      : Cost function parameter related to other elections
% vk    ((sizeXk-1) x4)  : Preference paramether.  Number of Xk x Number of Cand
% FAlph (2x1)        : Distribution of apha (beta dist), 2x1
% Sig_xsi (1x1)      : variance of Xsi
% DeltaO  (1x1)       : Increase of eligible voters for open election
% DeltaMO (1X1)       : Increase of eligible voters for modified opene elec.

% setting parameters
param(2,1)=-1;
C0       = 0;
Cx       = param(1:3,1);
Cz       = -1*abs(param(4:6,1));
FAlph(1,1)=abs(param(7,1));
FAlph(2,1)=abs(param(8,1));
param(7:8,1)= FAlph(1:2,1);
Sig_xsi  = max(0.5,abs(param(9,1)));
param(9,1) = Sig_xsi;
DeltaO  = 0.6891;
DeltaMO = 0.5366;
vk       = param(12:23,1);
composite= param(75:149,1);
Tij      = abs(param(150:260,1));
% new parameters
rho_eta = param(261,1);
rho_chi = param(262:265,1);
mu_chi = param(266:269,1);
chi = param(270:273,1)


% making log likelihood

# simulated values
N_T_ij=[0;1;3;6];        % The number of T_ij that we need when the number of candidates are 1,2,3,4.
Alpha    = icdf('beta',SimAlp,FAlph(1,1),FAlph(2,1));
Xsi      = Sig_xsi*SimXsi;

% calculating utilities
% Base Utility for Sincere, [# of categories X # of candidates]
VSin     = X(:,2:end)*[vk(1:3),vk(4:6),vk(7:9),vk(10:12)]-C0-X(:,2:end)*Cx(1:3)*ones(1,4);
% Base Utility for Strategic, [# of categories X # of candidates]
VStr     = X(:,2:end)*[vk(1:3),vk(4:6),vk(7:9),vk(10:12)];

%Eligible voters accounting for Open and Modified Open
DELTA    = DeltaO*Open+DeltaMO*MOpen;
RTOT     = RDemHat.*(1+DELTA)-VOther;

% store signals in advance
signals = randn(4,T) * sqrt(1/rho_eta)
for i = 1:4
  signals(i,:) = signals(i,:) + chi(i,1)

loglik_s=zeros(size(Cand,1),1);

% Candをdateについての昇順に変える
sort(Cand,13)

for S=1:size(Cand,1)
    if Cand(S,15)-Cand(S,14)<21 || S==30 || S==34 %Excluding Utah, Wisconsin, and small
    else
        N_candS = sum(Cand(S,1:4));
        M = Cand(S,15)-Cand(S,14)+1;      % Number of municipalities in State S
        if Cand(S,11)==0
            T_s=zeros(N_T_ij(N_candS,1),1);
            COMPOSITE=zeros(N_candS,1);
        else
            T_s=Tij(Cand(S,11):Cand(S,12),1);
            COMPOSITE=composite(Cand(S,16):Cand(S,17),1);
        end

        VSTR_s=[];
        Dropped_s=find(Cand(S,1:4)==0);            %Index of candidates withdrawn.
        Candidate_s = find(Cand(S,1:4)==1);        % index of cadidate
        Senate_s=DATA(Cand(S,14),27)*Cz(1,1);      %Cost from Senate elections
        Governer_s=DATA(Cand(S,14),29)*Cz(3,1);    %Cost from GOvernor elections
        dFX_s = dFX(Cand(S,14):Cand(S,15),:);
        date = Cand(S,13)                          % election date (1 ~ 14)

        temp=Cand(S,1:4);
        PatternCand_d=Cand(S,24:27);
        temp(:,(Cand(S,24:27)==0))=[];
        PatternCand_d(:,(Cand(S,24:27)==0))=[];

        Alpha_s=Alpha(Cand(S,14):Cand(S,15),:,:);
        Xsi_s=Xsi(Cand(S,14):Cand(S,15),:,:);

        VSin_s=VSin;
        VSin_s(:,Dropped_s)=[];
        VStr_s=VStr;
        VStr_s(:,Dropped_s)=[];

        % make XiOmg_s whose size is [1, num of candidates in the state]
        % take necessary parameters
        rho_chi_s = rho_chi(Candidate_s, 1)
        mu_chi_s = mu_chi(Candidate_s, 1)
        chi_s = chi_s(Candidate_s, 1)
        % calculating XiOmg_s
        signals_s = signals(:, date)
        cum_signals = cumsum(signals_s, 2)
        upper = rho_chi_s + mu_chi_s + rho_eta * cum_signals(Candidate_s,1)
        under = rho_chi_s + date * rho_eta
        XiOmg_s = upper ./ under

        % XiOmg_sの後ろのダッシュをひとつ消してください
        VSin_s=VSin_s+ones(size(VSin_s,1),1)*XiOmg_s.''-Senate_s-Governer_s;

        if Cand(S,9)==1 && Cand(S,10)==1
            Composite=ones(size(X,1),1)*COMPOSITE';
            VSTR_s(:,1)=T_s(1,1)*(VStr_s(:,1)-VStr_s(:,2))+T_s(2,1)*(VStr_s(:,1)-VStr_s(:,3))...
                +T_s(3,1)*(VStr_s(:,1)-VStr_s(:,4))+Composite(:,1);
            VSTR_s(:,2)=T_s(1,1)*(VStr_s(:,2)-VStr_s(:,1))+T_s(4,1)*(VStr_s(:,2)-VStr_s(:,3))...
                +T_s(5,1)*(VStr_s(:,2)-VStr_s(:,4))+Composite(:,2);
            VSTR_s(:,3)=T_s(2,1)*(VStr_s(:,3)-VStr_s(:,1))+T_s(4,1)*(VStr_s(:,3)-VStr_s(:,2))...
                +T_s(6,1)*(VStr_s(:,3)-VStr_s(:,4))+Composite(:,3);
            VSTR_s(:,4)=T_s(3,1)*(VStr_s(:,4)-VStr_s(:,1))+T_s(5,1)*(VStr_s(:,4)-VStr_s(:,2))...
                +T_s(6,1)*(VStr_s(:,4)-VStr_s(:,3))+Composite(:,4);

        elseif (Cand(S,9)==2 || Cand(S,9)==3) && Cand(S,10)==1
            Composite=ones(size(X,1),1)*COMPOSITE';
            VSTR_s(:,1)=T_s(1,1)*(VStr_s(:,1)-VStr_s(:,2))+T_s(2,1)*(VStr_s(:,1)-VStr_s(:,3))...
                +Composite(:,1);
            VSTR_s(:,2)=T_s(1,1)*(VStr_s(:,2)-VStr_s(:,1))+T_s(3,1)*(VStr_s(:,2)-VStr_s(:,3))...
                +Composite(:,2);
            VSTR_s(:,3)=T_s(2,1)*(VStr_s(:,3)-VStr_s(:,1))+T_s(3,1)*(VStr_s(:,3)-VStr_s(:,2))...
                +Composite(:,3);


        elseif Cand(S,10)==0  %after super tuesday
            % VSTR_s=VSin_s+C0+X(:,2:end)*Cx(1:3)*ones(1,N_candS) ?
            VSTR_s=VSin_s+2*(Senate_s+Governer_s);
        end

        % Utiltiy of Strategic with no house elections
        VSTR_s  = VSTR_s-C0-X(:,2:end)*Cx(1:3)*ones(1,N_candS)-Senate_s-Governer_s;

        VSTR_ss = zeros(N_sim,N_candS);
        VSIN_ss = zeros(N_sim,N_candS);

        % eligible voters
        RTot_s  = RTOT(Cand(S,14):Cand(S,15),:);
        RTot_s = max(RTot_s,sum(Votes(Cand(S,14):Cand(S,15),:),2));
        Votes_s = Votes(Cand(S,14):Cand(S,15),:)./(RTot_s*ones(1,4)); %vote share data % ここら辺はわかる
        Votes_s(:,Dropped_s)=[];
        loglik_m = zeros(M,1);

        if Cand(S,10)==1

            for m=1:M
                VSin_s=VSin_s-Cz(2,1)*DATA(Cand(S,14)+m-1,28);
                VSTR_s=VSTR_s-Cz(2,1)*DATA(Cand(S,14)+m-1,28);

                for sim=1:N_sim

                    nakami=max(min(VSin_s+ones(N_dFX,1)*Xsi_s(m,1:N_candS,sim),200),-200);
                    nakami2=max(min(VSTR_s+ones(N_dFX,1)*Xsi_s(m,1:N_candS,sim),200),-200);

                    eVSIN_ss=exp(nakami)./...
                        (1+sum(exp(nakami),2)*ones(1,N_candS));
                    eVSTR_ss=exp(nakami2)./...
                        (1+sum(exp(nakami2),2)*ones(1,N_candS));

                    VSTR_ss(sim,:)=dFX_s(m,:)*eVSTR_ss;
                    VSIN_ss(sim,:)=dFX_s(m,:)*eVSIN_ss;
                end

                Alp_ss = squeeze(Alpha_s(m,:,:)).'';
                VSHARE=VSTR_ss.*(Alp_ss*ones(1,N_candS))+VSIN_ss.*(1-Alp_ss*ones(1,N_candS));
                loglik_m(m,1)=log(sum(prod(pdf('norm',(ones(N_sim,1)*Votes_s(m,:)-VSHARE)/bandwidth,0,1),2),1)/N_sim);
            end
            loglik_s(S,1)=sum(loglik_m);

        elseif Cand(S,10)==0

            for m=1:M
                VSin_s=VSin_s-Cz(2,1)*DATA(Cand(S,14)+m-1,28);
                AST   =1./(1+exp(C0+X(:,2:end)*Cx(1:3)*ones(1,N_candS)...
                    +Cz(2,1)*DATA(Cand(S,14)+m-1,28)+Senate_s+Governer_s));
                %AST: Turnout of strategic voters after super tuesday

                for sim=1:N_sim

                    nakami=max(min(VSin_s+ones(N_dFX,1)*Xsi_s(m,1:N_candS,sim),200),-200);
                    eVSIN_ss=exp(nakami)./...
                        (1+sum(exp(nakami),2)*ones(1,N_candS));
                    VSIN_ss(sim,:)=dFX_s(m,:)*eVSIN_ss;
                end
                VSHARE=VSIN_ss;
                loglik_m(m,1)=log(sum(prod(pdf('norm',(ones(N_sim,1)*Votes_s(m,:)-VSHARE)/bandwidth,0,1),2),1)/N_sim);
            end
            loglik_s(S,1)=sum(loglik_m);
        end
    end
    if S==10
        S=S
    end
end

loglik=sum(loglik_s);
loglik=-loglik

n_eval=n_eval+1;
paramet=[loglik;n_eval;param];

file0=char(strcat({'parameter_'},{ii},{'_'},{jj},{'_'},{kk},{'.txt'}));
file1=char(strcat({'current_'},{ii},{'_'},{jj},{'_'},{kk},{'.txt'}));
current=paramet;
save(file1,'current','-ascii');
parameter=load(file0,'-ascii');
if bestpara(1,1)>paramet(1,1);
    bestpara=paramet;
    save(file0,'bestpara','-ascii');
end

if isnan(loglik)==1
    NANDERU=NaN
end
