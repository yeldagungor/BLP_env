    #T: # of markets; J: # of inside goods; 
    #K: # of measured product characteristics; L: measured marginal cost shifters

    #Measured Characteristics
    μVec_1 = ones(K + L); # Set mean 
    #ΣMat_1 = Matrix{Float64}(Diagonal(fill(0.9, (K + L, K + L))))+0.1*ones(K + L , K + L); #Set covariance matrix
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
    mc = exp.( 0.5 .+ charProduct_1[3,:,:] + 1.5 * charProduct_1[4,:,:] + charProduct_2[2,:,:] );

    #Random coefficients depend on integration methods:

    #part a
    include("numMethod_a.jl");
    S = 1000000; #1 million simulations
    priceVec = zeros(J,T); #Initialize
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
            @NLobjective(model, Min, 1)

            @NLconstraints(model, begin FOC1(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10) == 0
            FOC2(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10) == 0
            FOC3(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10) == 0
            FOC4(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10) == 0
            FOC5(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10) == 0
            FOC6(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10) == 0
            FOC7(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10) == 0
            FOC8(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10) == 0
            FOC9(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10) == 0
            FOC10(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10) == 0 end)
            
         
            p = [p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10];
            set_silent(model)
            JuMP.optimize!(model)
            @show JuMP.termination_status(model)
            priceVec[:,t] = JuMP.value.(p)
            
        end
    end
    


    #part b


    #part c

    #part d


    #part e

    #part f
