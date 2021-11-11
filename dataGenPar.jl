    #T: # of markets; J: # of inside goods; 
    #K: # of measured product characteristics; L: measured marginal cost shifters

    #Measured Characteristics
    μVec_1 = ones(K + L); # Set mean 
    ΣMat_1 = [1 0.1 0.1 0.1; 0.1 1 0.1 0.1; 0.1 0.1 0.25 0.1; 0.1 0.1 0.1 0.25]
    tempCharProduct_1 = MvNormal(μVec_1 , ΣMat_1);
    #Draw from MVnormal
    #First K rows belong to product char(meas); Last L rows belong to marginal cost shifters(meas) 

    #Columns belong to the product numbers; 3rd dimension is the market
    charProduct_1 = reshape(rand(tempCharProduct_1, T*J), (K + L , J , T));  

    #Unmeasured Characteristics
    μVec_2 = zeros(2); # Set mean 
    ΣMat_2 = [0.5 0.1;0.1 0.25]; #Set covariance matrix
    
    tempCharProduct_2 = MvNormal(μVec_2 , ΣMat_2);


    #Draw from MVnormal
    #First row belong to demand shifter(unmeas); Last row belong to marginal cost shifters(unmeas) 
    #Columns belong to the product numbers; 3rd dimension is the market
    charProduct_2 = reshape(rand(tempCharProduct_2, T*J), (2 , J , T));  

    #Marginal Costs
    #mc = exp.( 0.5 .+ charProduct_1[3,:,:] + 1.5 * charProduct_1[4,:,:] + charProduct_2[2,:,:] );
    mc = exp.( 0.5 .+ charProduct_1[3,:,:] + 1.5 * charProduct_1[4,:,:] + charProduct_2[2,:,:] );
    #mc = 0.5 .+ charProduct_1[3,:,:].^2 + 1.5 * charProduct_1[4,:,:].^2 + charProduct_2[2,:,:].^2 ;
    include("numMethod_a.jl");
    include("marketShare.jl");
    

    #part 2.a
    S = 1000000; #1 million simulations
    priceVec1 = zeros(J,T); #Initialize
    mShare = zeros(J+1,T)
    #beta draws
    d_0 = Normal(-2,1); β_0 = rand(d_0,(S,T)); 
    d_1 = Normal(1,0.5); β_1 = rand(d_1,(S,T));
    d_2 = Normal(1,0.5); β_2 = rand(d_2,(S,T));
    d_3 = Normal(1,0.5); β_3 = rand(d_3,(S,T));


    @time begin
        for k in 1:T
            global t = k
            initP = rand(1:10,10);

            model = Model(Ipopt.Optimizer);
            set_optimizer_attribute(model, "max_iter", 500)
            set_optimizer_attribute(model, "max_cpu_time",  3600.0)
            register(model,:FOC1, 10, FOC1; autodiff=true); #model object, function object, number of inputs, function
            register(model,:FOC2, 10, FOC2; autodiff=true);
            register(model,:FOC3, 10, FOC3; autodiff=true);
            register(model,:FOC4, 10, FOC4; autodiff=true);
            register(model,:FOC5, 10, FOC5; autodiff=true);
            register(model,:FOC6, 10, FOC6; autodiff=true);
            register(model,:FOC7, 10, FOC7; autodiff=true);
            register(model,:FOC8, 10, FOC8; autodiff=true);
            register(model,:FOC9, 10, FOC9; autodiff=true);
            register(model,:FOC10, 10, FOC10; autodiff=true); 

            @variable(model, mc[1,t] <= p_1 <=3000, start = mc[1,t] );
            @variable(model, mc[2,t] <= p_2 <=3000, start = mc[2,t] );
            @variable(model, mc[3,t] <= p_3 <=3000, start = mc[3,t] );
            @variable(model, mc[4,t] <= p_4 <=3000, start = mc[4,t] );
            @variable(model, mc[5,t] <= p_5 <=3000, start = mc[5,t] );
            @variable(model, mc[6,t] <= p_6 <=3000, start = mc[6,t] );
            @variable(model, mc[7,t] <= p_7 <=3000, start = mc[7,t] );
            @variable(model, mc[8,t] <= p_8 <=3000, start = mc[8,t] );
            @variable(model, mc[9,t] <= p_9 <=3000, start = mc[9,t] );
            @variable(model, mc[10,t] <= p_10 <=3000, start = mc[10,t] );
            @NLobjective(model, Min, 0)

            @NLconstraint(model, FOC1(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10) == 0)
           
            
            
            add_NL_constraint(model, :(FOC2($(p_1), $(p_2), $(p_3), $(p_4), $(p_5),$(p_6), $(p_7), $(p_8), $(p_9),$(p_10))== 0))
            add_NL_constraint(model, :(FOC3($(p_1), $(p_2), $(p_3), $(p_4), $(p_5),$(p_6), $(p_7), $(p_8), $(p_9),$(p_10))== 0))
            add_NL_constraint(model, :(FOC4($(p_1), $(p_2), $(p_3), $(p_4), $(p_5),$(p_6), $(p_7), $(p_8), $(p_9),$(p_10))== 0))
            add_NL_constraint(model, :(FOC5($(p_1), $(p_2), $(p_3), $(p_4), $(p_5),$(p_6), $(p_7), $(p_8), $(p_9),$(p_10))== 0))
            add_NL_constraint(model, :(FOC6($(p_1), $(p_2), $(p_3), $(p_4), $(p_5),$(p_6), $(p_7), $(p_8), $(p_9),$(p_10))== 0))
            add_NL_constraint(model, :(FOC7($(p_1), $(p_2), $(p_3), $(p_4), $(p_5),$(p_6), $(p_7), $(p_8), $(p_9),$(p_10))== 0))
            add_NL_constraint(model, :(FOC8($(p_1), $(p_2), $(p_3), $(p_4), $(p_5),$(p_6), $(p_7), $(p_8), $(p_9),$(p_10))== 0))
            add_NL_constraint(model, :(FOC9($(p_1), $(p_2), $(p_3), $(p_4), $(p_5),$(p_6), $(p_7), $(p_8), $(p_9),$(p_10))== 0))
            add_NL_constraint(model, :(FOC10($(p_1), $(p_2), $(p_3), $(p_4), $(p_5),$(p_6), $(p_7), $(p_8), $(p_9),$(p_10))== 0)) 

            p = [p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10];
            set_silent(model)
            JuMP.optimize!(model)
            @show JuMP.termination_status(model)
            
            priceVec1[:,t] = JuMP.value.(p)
            mShare[:,t]      = marketShare_1(t,S,priceVec1[:,t])[1];
        end
    end
    
    tempdf1 = DataFrame(P = ["P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"],
    T1 = priceVec1[:,1],T2 = priceVec1[:,2],T3 = priceVec1[:,3],
    T4 = priceVec1[:,4],T5 = priceVec1[:,5],T6 = priceVec1[:,6],
    T7 = priceVec1[:,7],T8 = priceVec1[:,8],T9 = priceVec1[:,9],T10 = priceVec1[:,10])
    tempdf2 = DataFrame(P = ["mean","med","std"],T1 = [mean(priceVec1[:,1]),median(priceVec1[:,1]),std(priceVec1[:,1])],T2 = [mean(priceVec1[:,2]),median(priceVec1[:,2]),std(priceVec1[:,2])],
    T3 = [mean(priceVec1[:,3]),median(priceVec1[:,3]),std(priceVec1[:,3])],T4 = [mean(priceVec1[:,4]),median(priceVec1[:,4]),std(priceVec1[:,4])],T5 = [mean(priceVec1[:,5]),median(priceVec1[:,5]),std(priceVec1[:,5])]
    ,T6 = [mean(priceVec1[:,6]),median(priceVec1[:,6]),std(priceVec1[:,6])],T7 = [mean(priceVec1[:,7]),median(priceVec1[:,7]),std(priceVec1[:,7])],T8 = [mean(priceVec1[:,8]),median(priceVec1[:,8]),std(priceVec1[:,8])],
    T9 = [mean(priceVec1[:,9]),median(priceVec1[:,9]),std(priceVec1[:,9])],T10 = [mean(priceVec1[:,10]),median(priceVec1[:,10]),std(priceVec1[:,10])])
    df1 = vcat(tempdf1,tempdf2)
        

    tempdf1 = DataFrame(P = ["P0","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"],
    T1 = mShare[:,1],T2 = mShare[:,2],T3 = mShare[:,3],
    T4 = mShare[:,4],T5 = mShare[:,5],T6 = mShare[:,6],
    T7 = mShare[:,7],T8 = mShare[:,8],T9 = mShare[:,9],T10 = mShare[:,10])
    tempdf2 = DataFrame(P = ["mean","med","std"],T1 = [mean(mShare[2:11,1]),median(mShare[2:11,1]),std(mShare[2:11,1])],T2 = [mean(mShare[2:11,2]),median(mShare[2:11,2]),std(mShare[2:11,2])],
    T3 = [mean(mShare[2:11,3]),median(mShare[2:11,3]),std(mShare[2:11,3])],T4 = [mean(mShare[2:11,4]),median(mShare[2:11,4]),std(mShare[2:11,4])],T5 = [mean(mShare[2:11,5]),median(mShare[2:11,5]),std(mShare[2:11,5])]
    ,T6 = [mean(mShare[2:11,6]),median(mShare[2:11,6]),std(mShare[2:11,6])],T7 = [mean(mShare[2:11,7]),median(mShare[2:11,7]),std(mShare[2:11,7])],T8 = [mean(mShare[2:11,8]),median(mShare[2:11,8]),std(mShare[2:11,8])],
    T9 = [mean(mShare[2:11,9]),median(mShare[2:11,9]),std(mShare[2:11,9])],T10 = [mean(mShare[2:11,10]),median(mShare[2:11,10]),std(mShare[2:11,10])])
    df2 = vcat(tempdf1,tempdf2)
        
  
