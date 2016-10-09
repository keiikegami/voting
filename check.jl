# define T as num of period
T = 14
bandwidth = 0.05
N_sim = 100
n_eval = 0
I = 4
N_cand = 4
# load iii.txt -ASCII;
# load jjj.txt -ASCII;

# dataのロード
data = readtable("data.csv", header = false)
data = convert(Array, data)

#make file
iii = 1 #kokokaeruiii;
jjj = 1 #kokokaerujjj;
kkk = 2 #kokokaerukkk;

# convert int to string
ii = "$iii";
jj = "$jjj";
kk = "$kkk";

# たぶんココいらない
file3 = *("parameter_", "$ii", "_", "$jj", "_", "$kk", ".txt")
file4 = *("inivalue_", "$ii", ".txt")

# これはinivalueをロードしてる
inivalue = readtable(*("inivalue_", "$ii", ".txt"), separator = '\t')

inivalue = inivalue[:, jjj]

# 直下で定義したparameterをfile3の名前で保存
parameter = [1000000000;0;inivalue]
writetable(*("parameter_", "$ii", "_", "$jj", "_", "$kk", ".txt"), DataFrame(X = parameter), separator = '\t')

#　ここから本題
D　=　data[:,45] .> 100 # Remove municipality with small populaion <100
# 行削除はどうするのか、そもそもいるのか
# data[D, :] = [];
data = data[D, :]

# Construct X : Regressor for individual characteristics
X = [1 1 1 1 1];
dFX = ones(size(data, 1))
Race = [1 0 0; 0 1 0; 0 0 1];
Educ = 0.1*[16;14;12;9];
Incm = log([20000;35000;72500;120000]);

# (white,black,otherasian+indian+other)
dFXRace = [data[:,121] data[:,124] data[:,122] + data[:,123] + data[:,125]];

dFXEduc = [data[:,126:128] (1- sum(data[:,126:128], 2))];  # (overba,underba,hs,other)
dFXIncm = [sum(data[:,129:132], 2) sum(data[:,133:136], 2) sum(data[:,137:140], 2) sum(data[:,141:144], 2)];

# dFXIncm=data(:,129:144); % (income 16 category:
# -10K,15K,20K,25K,30K,35K,40,45,50,60,75,100,125,150,200,200+)
# ここは精査する必要がありそう。
for r in 1:3
    for e in 1:4
        for i in 1:4
            X = [X; reshape([Race[1,1:3] ; Educ[1] ; Incm[1]], 1, 5)];
            dFX = [dFX dFXRace[:,r].*dFXEduc[:,e].*dFXIncm[:,i]];
        end
    end
end

del = trues(size(X, 1))
del[1] = false
X = X[del, :]

N_dFX = size(dFX,2);
N_muni = size(dFX,1);      # (number of municipalities)
# data3,4,5,6はまさかのstring型
Votes = [max(0,data[:,3]) max(0,data[:,4]) max(0,data[:,5]) max(0,data[:,6])]; # (votes:  clark, dean, edwards, kerry)municipalityごとに、いない奴は0になるようにしてる
VTotal = data[:,117];      # (total number of votes)
RDemHat = data[:,148];     # (registred number of democrats predicted)
PopTot = data[:,45];       # (total population)
Open = data[:,118];        # open
MOpen = data[:,119];       # modified open
VOther = VTotal - sum(Votes,2);  #Votes of Penna & Sharpton etc.
RegDem = RDemHat./PopTot;  # (fraction of registered democrats in population)

# make Cand
Cand = zeros(35,17);       

                         # Candidates whose name is on ballot (column 1;4)
                         # Column 5 corresponds to ichiban maeno parameter
                         # noichi (for T_ij), Col6 to ichiban saigono parameter.
                         # Same for Columns 7 and 8 (for E[qi|omega]).
                         # Column 9: Who's in the race.
                         # Column 10: 1 if Before (and including) super tuesday
                         # Column 11:　？
                         # Column 12:　？
                         # Column 13: Date t (1,2,...14). 6 correspond to
                         # super Tues.
                         # Column 14: Municipality Number Start
                         # Column 15: Municipality Number End;
                         # 下のE[qi|omega_piv]-E[qi|omega_piv]って0では？
                         # Column 16: Atamadashi for
                         # E[qi|omega_piv]-E[qi|omega_piv]
                         # Column 17: Owaridashi for
                         # E[qi|omega_piv]-E[qi|omega_piv]
                         # the remains is explained below

j = 0
N_T_ij = [0;1;3;6];        # The number of T_ij that we need when the number of candidates are 1,2,3,4.
for i in 1:size(data, 1)
    if j == data[i,149] - 1
        Cand[j+1,1:4] = squeeze(ones(1,4), 1)+min(data[i,3:6], squeeze(zeros(1,4), 1));
        if j == 0
            Cand[j+1,5] = 0;
            Cand[j+1,6] = 0;
            Cand[j+1,11] = 0;
        else
            Cand[j+1,5] = Cand[j,6] + (sum(Cand[j+1,1:4])>1);
            Cand[j+1,6] = Cand[j+1,5] + N_T_ij[convert(Int64, sum(Cand[j+1,1:4])),1] - 1+(sum(Cand[j+1,1:4])==1);
        end

        # 1 if before or on super Tues.

        Cand[j+1,10] = (data[i,120] - 92<=0);
        Cand[j+1,7] = sum(sum(Cand[1:j+1,1:4]))-sum(Cand[j+1,1:4])+1;
        Cand[j+1,8] = sum(sum(Cand[1:j+1,1:4]));

        Cand[j+1,17] =  sum(sum(Cand[1:j+1,1:4].*(Cand[1:j+1,10]*ones(1,4))));
        Cand[j+1,16] = Cand[j+1,17]-sum(Cand[j+1,1:4]*Cand[j+1,10])+1;
        Cand[j+1,16:17] = Cand[j+1,16:17]*Cand[j+1,10];

        # Baai-wake: [Clark,Dean,Edawards]
        # Type 0: [000] Type 1:[111] Type2:[101] Type3: [011] Type4: [010]
        if Cand[j+1,1:3] == zeros(1,3)
            Cand[j+1,9] = 0;
        elseif Cand[j+1,1:3] == ones(1,3)
            Cand[j+1,9] = 1;
        elseif Cand[j+1,1:3] == [1,0,1]
            Cand[j+1,9] = 2;
        elseif Cand[j+1,1:3] == [0,1,1]
            Cand[j+1,9] = 3;
        else
            Cand[j+1,9] = 4;
        end

        #Atamadashi for params
        if j > 0
            Cand[j+1,11] = sum((Cand[1:j,6] - Cand[1:j,5]+1).*Cand[1:j,10],1)+Cand[j+1,10]*(sum(Cand[j+1,1:4]) > 1);
            Cand[j+1,12] = sum((Cand[1:j+1,6] - Cand[1:j+1,5]+1).*Cand[1:j+1,10],1);
        end
        
        Cand[j+1,13] = data[i,120];
        Cand[j+1,14] = i;
         j = j+1;
        
    end


    Cand[:,11] = Cand[:,10].*Cand[:,11];      # Atamadashi ignoring post-super
    Cand[:,12] = Cand[:,10].*Cand[:,12];      # Tues states
    Cand[1:end-1,15] = Cand[2:end,14]-1;
    Cand[end,15] = i;

end

iij = unique(Cand[:,13]);
for i in 1:size(unique(Cand[:,13]),1)
    for j in 1:size(Cand[:,13],1)
        if Cand[j,13] == iij[i]
            Cand[j,13] = i;
        end
    end
end

NNCan[1,1] = 4;
A_Exi[1,:] = [1 4];
PatternCandall[1,:] = [1 1 1 1];
for S in 2:max(Cand[:,13], 0, 1)
    # of candidate on date S
    Temp = Cand;
    # 行削除はどうする
    # Temp[(Temp[:,13] .!= S), :]=[];
    Temp = Temp[[Temp[:,13] .== S], :]
    NNCan[S,1] = sum((sum(Temp[:,1:4],1) .> 0),2);
    PatternCandall[S,1:4] = (sum(Temp[:,1:4],1) > 0);
    # Atamadashi for ExiOmega
    A_Exi[S,:] = [A_Exi[S-1,2] + 1, A_Exi[S-1,2] + NNCan[S,1]];
end

for S = 1:size(Cand,1)
    SS = Cand[S,13];
    Cand[S,22:23] = [A_Exi[SS,:]];
    Cand[S,24:27] = PatternCandall[SS,:];
end


# parameter description
# Xi|Omg (118x1)     : For each state, we have at most 4 values (total 118)
# Xi|OmgPiv (56x1)   : For each state, we have at most 3 values (total  56)
# Tij (105x1)        : For each state, we have at most 6 values (total 105)
#  Set to zero if after Super Tuesday.
# Cx    (sizeXk x1)  : Cost function parameter.  Number of Xk
# Cz    (3 x 1)      : Cost function parameter related to other elections
# vk    (sizeXk x4)  : Preference paramether.  Number of Xk x Number of Cand
# FAlph (2x1)        : Distribution of alpha (beta dist), 2x1
# Sig_xsi (1x1)      : variance of Xsi
# DeltaO  (1x1)       : Increase of eligible voters for open election
# DeltaMO (1X1)       : Increase of eligible voters for modified opene elec.

# load inivalue.txt inivalue -ASCII
# load indicator.txt indicator -ASCII;
# load inival.txt inival -ASCII

DATA = data;
srand(10);

# add XiOmg
# load the simulated XiOmg data
# reshape the data to [T, N_cand, N_sim]
SimXsi = randn(N_muni, N_cand, N_sim);
SimAlp = rand(N_muni, N_sim);

# simulate for SimXiOmg
# initial values are maked by python
# learning_params is [13, 100]
# maybe use just the first column
learning_params = readtable('learning_params.txt', separator = '\t')
# append learing_params(:,1) to the last of inivalue
inivalue = [inivalue; learning_params[:,1]]


# delete OPTIONS1
# replace x0 by inivalue
bestpara = [1000000;0;inivalue];
theta, likli = optimize(new_loglike, inivalue, BFGS())
