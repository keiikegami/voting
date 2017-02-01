# preparation

using Optim
using StatsFuns
using DataFrames
using Gadfly
using PyPlot

import StatsFuns.betainvcdf
betainvcdf(alpha::Number, beta::Number, x::Array) = reshape([betainvcdf(alpha, beta, i) for i in reshape(x, 1, size(x, 1)*size(x, 2)) ], size(x, 1), size(x, 2))
import Base.max
max(number::Real, comparison::Any) = [max(number, parse(Int, i)) for i in comparison] 
import StatsFuns.normpdf
normpdf(array::Array{Float64, 2}) = reshape([normpdf(i) for i in reshape(array, 1, size(array, 1)*size(array, 2))], size(array, 1), size(array, 2))

include("bayes.jl");


# type define
type simulation
    param::Array{Float64, 1}
end


# methods are defined as follows

# method 1 : simulate vote rates
function simulate(m::simulation)
    simu_vote = Array(Float64, size(DATA, 1), 4)

    # パラメータセット
    FAlph = Array(Float64, 2, 1)
    param = m.param
    FAlph = Array(Float64, 2, 1)
    # setting parameters
    param[2] = -1
    C0 = 0
    # Raceの変更に伴い、Cxも4つの要素にする必要がある。
    q = 4
    Cx = param[1:q]
    Cz = -abs(param[q+1:q+3])
    FAlph[1,1] = abs(param[q+4])
    FAlph[2,1] = abs(param[q+5])
    param[q+4:q+5] = FAlph[1:2,1]
    Sig_xsi  = max(0.5, abs(param[q+6]))
    param[q+6] = Sig_xsi
    DeltaO  = 0.6891
    DeltaMO = 0.5366
    # Raceをいじった関係でこのvkを4*4 = 16この要素にしなくちゃいけない
    vk = param[12:27]
    composite = param[75:149]
    Tij = abs(param[150:260])

    # new parameters
    # rho_eta = abs(param[261])
    rho_eta = 1
    rho_chi = param[262:265]
    mu_chi = param[266:269]
    chi = param[270:273]


    N_T_ij = [0;1;3;6]
    # 上で拡張した関数を使用
    Alpha = betainvcdf(FAlph[1,1],FAlph[2,1], SimAlp)
    Xsi = Sig_xsi*SimXsi

    # calculating utilities
    # Base Utility for Sincere, [# of categories X # of candidates]
    # VSinは48×4
    VSin = X[:,2:end]*[vk[1:4] vk[5:8] vk[9:12] vk[13:16]]-C0-X[:, 2:end]*Cx*ones(1,4)
    # Base Utility for Strategic, [# of categories X # of candidates]
    # VStrも48×4
    VStr = X[:, 2:end]*[vk[1:4] vk[5:8] vk[9:12] vk[13:16]]

    #Eligible voters accounting for Open and Modified Open
    DELTA = DeltaO*Open+DeltaMO*MOpen
    RTOT = RDemHat.*(1+DELTA)-VOther

    # store signals in advance
    signals = randn(4,T) * sqrt(1/rho_eta)

    for i in 1:4
        signals[i, :] = signals[i, :] + chi[i, 1]
    end

    for S in 1:size(Cand,1)
        #println(S)

        if Cand[S, 15] - Cand[S, 14] < 21 || S == 30 || S == 34 #Excluding Utah, Wisconsin, and small
        else

            N_candS = sum(Cand[S,1:4] .!= 0)
            M = Cand[S, 15] - Cand[S, 14] + 1    # Number of municipalities in State S

            # 州ごとにmunicipality別の得票率を収める箱を作る
            # dropしてる人のところに999とか入れて、municipality数*4のArrayにぶちこむ
            # 本来のループではそれを盾に結合して全部の州についてまとまってるArrayを作る。
            simu_votes = Array(Float64, M, 4)

            if Cand[S, 11] == 0
                T_s = zeros(N_T_ij[N_candS,1],1)
                COMPOSITE = zeros(N_candS,1)
            else
                T_s = Tij[Cand[S, 11]:Cand[S, 12],1]
                COMPOSITE = composite[Cand[S, 16]:Cand[S, 17],1]
            end

            # VSTR_s = []
            Dropped_s = find(Cand[S,1:4] .== 0)            # index of candidates withdrawn.
            Rem = Cand[S,1:4] .!= 0
            Candidate_s = find(Cand[S,1:4] .== 1)        # index of cadidate
            Senate_s = DATA[Cand[S, 14], 27]*Cz[1,1]     # Cost from Senate elections
            Governer_s = DATA[Cand[S, 14], 29]*Cz[3,1]    # Cost from GOvernor elections
            dFX_s = dFX[Cand[S, 14]:Cand[S, 15], :]
            date = Cand[S, 13]                          # election date (1 ~ 14)

            # というか列の削除をする必要があるのか
            temp = Cand[S,1:4]
            PatternCand_d = Cand[S, 24:27]
            # temp[:, find(Cand[S,24:27] .== 0)] = []
            # PatternCand_d[:, find(Cand[S,24:27] .== 0)] = []
            temp = temp[Cand[S,24:27] .!= 0.0]
            PatternCand_d = PatternCand_d[Cand[S,24:27] .!= 0.0]

            Alpha_s = reshape(Alpha[Cand[S, 14]:Cand[S, 15], :, :], length(Cand[S, 14]:Cand[S, 15]), 100)
            Xsi_s = Xsi[Cand[S, 14]:Cand[S, 15], :, :]

            # そもそも削除する必要があるか
            VSin_s = VSin
            # VSin_s[:, Dropped_s] = []
            VSin_s = VSin_s[:, Rem]
            VStr_s = VStr
            # println(size(VStr_s, 2))
            # VStr_s[:, Dropped_s] = []
            VStr_s = VStr_s[:, Rem]

            # make XiOmg_s whose size is [1, num of candidates in the state]
            # take necessary parameters
            rho_chi_s = rho_chi[Candidate_s, 1]
            mu_chi_s = mu_chi[Candidate_s, 1]
            chi_s = chi[Candidate_s, 1]
            # calculating XiOmg_s
            # signals_s = signals[:, date]

            # cumsum は1次元配列にしか使えないので注意ダメだったらsqueeze
            cum_signals = cumsum(signals, 2)
            upper = rho_chi_s + mu_chi_s + rho_eta * cum_signals[Candidate_s, date]
            under = rho_chi_s + date * rho_eta
            XiOmg_s = upper ./ under

            VSin_s = VSin_s + ones(size(VSin_s,1), 1)*XiOmg_s' -Senate_s-Governer_s

            # col9は0~4、4は？

            if Cand[S, 9] == 1 && Cand[S, 10] == 1
                VSTR_s = ones(48, 4)
                Composite = ones(size(X,1),1)*COMPOSITE'
                VSTR_s[:,1] = T_s[1,1]*(VStr_s[:,1]-VStr_s[:,2])+T_s[2,1]*(VStr_s[:,1]-VStr_s[:,3])
                +T_s[3,1]*(VStr_s[:,1]-VStr_s[:,4])+Composite[:,1]
                VSTR_s[:,2] = T_s[1,1]*(VStr_s[:,2]-VStr_s[:,1])+T_s[4,1]*(VStr_s[:,2]-VStr_s[:,3])
                +T_s[5,1]*(VStr_s[:,2]-VStr_s[:,4])+Composite[:,2]
                VSTR_s[:,3] = T_s[2,1]*(VStr_s[:,3]-VStr_s[:,1])+T_s[4,1]*(VStr_s[:,3]-VStr_s[:,2])
                +T_s[6,1]*(VStr_s[:,3]-VStr_s[:,4])+Composite[:,3]
                VSTR_s[:,4] = T_s[3,1]*(VStr_s[:,4]-VStr_s[:,1])+T_s[5,1]*(VStr_s[:,4]-VStr_s[:,2])
                +T_s[6,1]*(VStr_s[:,4]-VStr_s[:,3])+Composite[:,4]

            elseif (Cand[S, 9] == 2 || Cand[S, 9] == 3 || Cand[S, 9] == 4) && Cand[S, 10] == 1
                VSTR_s = ones(48, 3)
                Composite = ones(size(X,1),1)*COMPOSITE'
                VSTR_s[:,1] = T_s[1,1]*(VStr_s[:,1]-VStr_s[:,2])+T_s[2,1]*(VStr_s[:,1]-VStr_s[:,3])
                +Composite[:,1]
                VSTR_s[:,2] = T_s[1,1]*(VStr_s[:,2]-VStr_s[:,1])+T_s[3,1]*(VStr_s[:,2]-VStr_s[:,3])
                +Composite[:,2]
                VSTR_s[:,3] = T_s[2,1]*(VStr_s[:,3]-VStr_s[:,1])+T_s[3,1]*(VStr_s[:,3]-VStr_s[:,2])
                +Composite[:,3]


            elseif Cand[S, 10] == 0  # after super tuesday
                VSTR_s=VSin_s+C0+X[:,2:end]*Cx*ones(1,N_candS)
                # VSTR_s = VSin_s + 2*(Senate_s+Governer_s)
            end

            # Utiltiy of Strategic with no house elections
            VSTR_s = VSTR_s - C0 - X[:,2:end]*Cx*ones(1,N_candS) - Senate_s - Governer_s

            # eligible voters
            RTot_s = RTOT[Cand[S, 14]:Cand[S, 15], :]
            RTot_s = max(RTot_s, sum(Votes[Cand[S, 14]:Cand[S, 15], :], 2))
            Votes_s = Votes[Cand[S, 14]:Cand[S, 15], :]./(RTot_s*ones(1,4)) #vote share data
            # 行削除
            # Votes_s[:, Dropped_s] = []
            Votes_s = Votes_s[:, Rem]
            loglik_m = zeros(M,1)

            VSTR_ss = Array(Float64, N_sim,N_candS)
            VSIN_ss = Array(Float64, N_sim,N_candS)

            # ここ以降が遅い
            # どっちも遅いが、beforeの方が平均的に遅い
            if Cand[S, 10] == 1

                A1 = Array(Float64, N_dFX, N_candS)
                A2 = Array(Float64, N_dFX, N_candS)
                B1 = Array(Float64, 1, N_candS)
                B2 = Array(Float64, 1, N_candS)
                cont1 = Array(Float64, N_dFX)
                cont2 = Array(Float64, N_dFX)
                #eVSIN_ss = Array(Float64, N_dFX, N_candS)
                #eVSTR_ss = Array(Float64, N_dFX, N_candS)

                for m in 1:M
                    VSin_s = VSin_s - Cz[2,1]*DATA[Cand[S, 14] + m - 1, 28]
                    VSTR_s = VSTR_s - Cz[2,1]*DATA[Cand[S, 14] + m - 1, 28]

                    # for文に書き換え
                    for sim in 1:N_sim

                        # nakami = max(min(VSin_s + ones(N_dFX,1)*reshape(Xsi_s[m,1:N_candS, sim], 1, N_candS), 200), -200)
                        # nakami2 = max(min(VSTR_s + ones(N_dFX,1)*reshape(Xsi_s[m,1:N_candS,sim], 1, N_candS), 200), -200)
                        for i in 1:N_candS
                            for j in 1:N_dFX
                                A1[j, i] = VSin_s[j, i] + Xsi_s[m, i, sim]
                                A2[j, i] = VSTR_s[j, i] + Xsi_s[m, i, sim]
                            end
                        end
                        nakami = max(min(A1, 200.0), -200.0)
                        nakami2 = max(min(A2, 200.0), -200.0)


                        eVSIN_ss = exp(nakami)./ (1+sum(exp(nakami),2)*ones(1,N_candS))
                        eVSTR_ss = exp(nakami2)./ (1+sum(exp(nakami2),2)*ones(1,N_candS))
                        #naka = sum(exp(nakami), 2)
                        #naka2 = sum(exp(nakami2), 2)
                        #for i in 1:N_candS
                          #  for j in 1:N_dFX
                            #    eVSIN_ss[j, i] = exp(nakami[j,i])/(1+ naka[j])
                               # eVSTR_ss[j, i] = exp(nakami2[j,i])/(1+ naka2[j])
                            #end
                        #end


                        # VSTR_ss[sim, :] = reshape(dFX_s[m, :],1, 48)*eVSTR_ss
                        # VSIN_ss[sim, :] = reshape(dFX_s[m, :], 1, 48)*eVSIN_ss
                        for i in 1:N_candS
                            for j in 1:N_dFX
                                cont1[j] = dFX_s[m, j]*eVSTR_ss[j, i]
                                cont2[j] = dFX_s[m, j]*eVSIN_ss[j, i]
                            end
                            VSTR_ss[sim,  i] = sum(cont1)
                            VSIN_ss[sim, i] = sum(cont2)
                        end

                    end

                    Alp_ss = Alpha_s[m, :]
                    VSHARE = VSTR_ss.*(Alp_ss*ones(1,N_candS)) + VSIN_ss.*(1-Alp_ss*ones(1,N_candS))
                    simu_votes[m, Candidate_s] = mean(VSHARE, 1)
                    simu_votes[m, Dropped_s] = 0
                end
                simu_vote[Cand[S, 14]:Cand[S, 15], :] = simu_votes

            elseif Cand[S, 10] == 0

                A1 = Array(Float64, N_dFX, N_candS)
                B2 = Array(Float64, 1, N_candS)
                cont2 = Array(Float64, N_dFX)
                #eVSIN_ss = Array(Float64, N_dFX, N_candS)

                for m in 1:M
                    VSin_s = VSin_s-Cz[2,1]*DATA[Cand[S, 14]+m-1, 28]
                    AST = 1./(1+exp(C0+X[:,2:end]*Cx*ones(1,N_candS) + Cz[2,1]*DATA[Cand[S, 14]+m-1, 28]+Senate_s+Governer_s))
                    # AST: Turnout of strategic voters after super tuesday
                    for sim in 1:N_sim

                        # nakami = max(min(VSin_s + ones(N_dFX,1)*reshape(Xsi_s[m,1:N_candS, sim], 1, N_candS), 200), -200)
                        for i in 1:N_candS
                            for j in 1:N_dFX
                                A1[j, i] = VSin_s[j, i] + Xsi_s[m, i, sim]
                            end
                        end
                        nakami = max(min(A1, 200.0), -200.0)


                        eVSIN_ss = exp(nakami)./ (1 + sum(exp(nakami),2)*ones(1,N_candS))
                        # naka = sum(exp(nakami), 2)
                        # for i in 1:N_candS
                           # for j in 1:N_dFX
                              #  eVSIN_ss[j, i] = exp(nakami[j,i])/(1+ naka[j])
                            #end
                        #end


                        # VSIN_ss[sim, :] = reshape(dFX_s[m, :], 1, 48) * eVSIN_ss
                        for i in 1:N_candS
                            for j in 1:N_dFX
                                cont2[j] = dFX_s[m, j]*eVSIN_ss[j, i]
                            end
                            VSIN_ss[sim, i] = sum(cont2)
                        end

                    end
                    VSHARE = VSIN_ss
                    simu_votes[m, Candidate_s] = mean(VSHARE, 1)
                    simu_votes[m, Dropped_s] = 0
                end
                simu_vote[Cand[S, 14]:Cand[S, 15], :] = simu_votes
            end
        end
    end
    return simu_vote
end

# method 2 : visualization type 1 (for each combination of cadidates)
# default setting allows you to draw the results of all states
# you can set t as state number (1 ~ 35)
function candcand(m::simulation, t = 0, Votes = false)
    if Votes == false
        Votes = simulate(m)
    end
    
    if t == 0
        plt = PyPlot
        names = ["clark", "dean", "edwards", "kerry"]
        DeltaO  = 0.6891
        DeltaMO = 0.5366
        combination = [(i, i+j) for i in 1:3 for j in 1:(4-i)]

        for S in 1:size(Cand,1)

            DELTA = DeltaO*Open+DeltaMO*MOpen
            RTOT = RDemHat.*(1+DELTA)-VOther
            RTot_s = RTOT[Cand[S, 14]:Cand[S, 15], :]
            RTot_s = max(RTot_s, sum(Votes[Cand[S, 14]:Cand[S, 15], :], 2))
            Votes_s = Votes[Cand[S, 14]:Cand[S, 15], :]./(RTot_s*ones(1,4))

            num_rows, num_cols = 2, 3
            fig, axes = subplots(num_rows, num_cols, figsize=(12, 8))
            axes = vec(axes)

            # cand1 vs cand2で、cand1が横軸、cand2が縦軸
            for (n,c) in enumerate(combination)
                cand1 = names[c[1]]
                cand2 = names[c[2]]
                ax = axes[n]
                ax[:scatter](Votes_s[:, c[1]], Votes_s[:, c[2]], s = 3)
                ax[:set_title]("$cand1 vs $cand2")
                #ax[:set_xticks]([0,0.25,0.5,0.75])
                #ax[:set_yticks]([0,0.25,0.5,0.75])
                savefig("state_number_$S")
            end
        end
    
    else
        plt = PyPlot
        names = ["clark", "dean", "edwards", "kerry"]
        DeltaO  = 0.6891
        DeltaMO = 0.5366
        combination = [(i, i+j) for i in 1:3 for j in 1:(4-i)]

        S = t
        DELTA = DeltaO*Open+DeltaMO*MOpen
        RTOT = RDemHat.*(1+DELTA)-VOther
        RTot_s = RTOT[Cand[S, 14]:Cand[S, 15], :]
        RTot_s = max(RTot_s, sum(Votes[Cand[S, 14]:Cand[S, 15], :], 2))
        Votes_s = Votes[Cand[S, 14]:Cand[S, 15], :]./(RTot_s*ones(1,4))

        num_rows, num_cols = 2, 3
        fig, axes = subplots(num_rows, num_cols, figsize=(12, 8))
        axes = vec(axes)

        # cand1 vs cand2で、cand1が横軸、cand2が縦軸
        for (n,c) in enumerate(combination)
            cand1 = names[c[1]]
            cand2 = names[c[2]]
            ax = axes[n]
            ax[:scatter](Votes_s[:, c[1]], Votes_s[:, c[2]], s = 3)
            ax[:set_title]("$cand1 vs $cand2")
            #ax[:set_xticks]([0,0.25,0.5,0.75])
            #ax[:set_yticks]([0,0.25,0.5,0.75])
            savefig("state_number_$S")
        end
    end
end

# method 3 : visualization type 2 (summary vote rate for each state)
function state_vote(m::simulation, Votes = false)
    
    DeltaO  = 0.6891
    DeltaMO = 0.5366
    DELTA = DeltaO*Open+DeltaMO*MOpen
    RTOT = RDemHat.*(1+DELTA)-VOther
    names = ["clark", "dean", "edwards", "kerry"]
    
    if Votes == false
        Votes = simulate(m)
    end
    
    shares = Array(Float64, size(Cand,1), 4)

    # municipalityごとにそのdemocratic人口で得票率をかけて、候補者ごとにその得票数を足し算。その後、その州におけるdemocratic人口の合計で割る。
    for S in 1:size(Cand, 1)
        shares[ S, :] = sum(Votes[Cand[S, 14]:Cand[S, 15], :] .* RTOT[Cand[S, 14]:Cand[S, 15], :], 1) ./ sum(RTOT[Cand[S, 14]:Cand[S, 15], :])
    end

    fig, ax = subplots()
    for i in 1:4
        cand = names[i]
        ax[:plot](shares[:, i], linewidth=2, alpha=0.6, label="$cand")
    end
    ax[:legend]()
    savefig("vote_share_states")
end

# method 4 : visualization type 3 (demographic scatter)
# you can choose demographic factor from ["race", "edu", "income"]
# default setting allows you to draw the result of all states
# you can set t as state number (1 ~ 35)
function demo(m::simulation, demogra::String, t = 0, Votes = false)
    if Votes == false
        Votes = simulate(m)
    end
    
    plt = PyPlot
    names = ["clark", "dean", "edwards", "kerry"]
    DeltaO  = 0.6891
    DeltaMO = 0.5366
    DELTA = DeltaO*Open+DeltaMO*MOpen
    RTOT = RDemHat.*(1+DELTA)-VOther

    voterate = sum(Votes, 2)./RTOT
    if demogra == "race"
        a = dFXEduc[:, 1]./sum(dFXEduc[:,2:4],2)
    elseif demogra == "edu"
        a = dFXRace[:, 1]./sum(dFXRace[:,2:3],2)
    elseif demogra == "income"
        a = dFXIncm[:, 1]
    end
    kerry = Votes[:, 4]./RTOT
    
    if t == 0
        num_rows, num_cols = 7, 5
        fig, axes = subplots(num_rows, num_cols, figsize=(12, 20))
        axes = vec(axes)
        
        for S in 1:size(Cand,1)

            voterate_s = voterate[Cand[S, 14]:Cand[S, 15], :]
            b = a[Cand[S, 14]:Cand[S, 15], :]
            kerry_s = kerry[Cand[S, 14]:Cand[S, 15], :]

            ax = axes[S]
            ax[:scatter](b, voterate_s, s = 50*kerry_s, alpha = 0.5)
            ax[:set_title]("state_number_$S")
            #ax[:set_xticks]([0,0.25,0.5,0.75])
            #ax[:set_yticks]([0,0.25,0.5,0.75])

        end
        # matplotlibそのまんま使えるのね
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        savefig("$demogra_plot")
        
    else
        S = t
        voterate_s = voterate[Cand[S, 14]:Cand[S, 15], :]
        b = a[Cand[S, 14]:Cand[S, 15], :]
        kerry_s = kerry[Cand[S, 14]:Cand[S, 15], :]

        plt.scatter(b, voterate_s, s = 50*kerry_s, alpha = 0.5)
        plt.title("state_number_$S")
        savefig("$demogra" *"_"* "$S"*"_plot")
    end
end

# method5 : make histogram of vote rate for each municipality
function rate_histo(m::simulation, Votes = false)
    if Votes == false
        Votes = simulate(m)
    end

    DeltaO  = 0.6891
    DeltaMO = 0.5366
    DELTA = DeltaO*Open+DeltaMO*MOpen
    RTOT = RDemHat.*(1+DELTA)-VOther
    RTOT = max(RTOT, sum(Votes, 2))
    Vote_rate = sum(Votes ./ RTOT, 2)
    Vote_rate = convert(DataFrame, Vote_rate)
    deleterows!(Vote_rate, find(isna(Vote_rate[:,1])))
    Gadfly.plot(Vote_rate,x = "x1",  Geom.histogram(bincount=100))
    
end