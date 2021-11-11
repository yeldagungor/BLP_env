function marketShare_1(t, S, price)
    β_0_t = β_0[:,t]; #size S
    β_1_t = β_1[:,t]; #size S
    β_2_t = β_2[:,t]; #Size S
    β_3_t = β_3[:,t]; #Size S
    ξ_t = charProduct_2[1,:,t];

    temp0 = zeros(S,J)
    temp1= ones(S,J)
    outContr = ones(S) #Outside option contribution 
    β_0_t_temp = β_0_t.*temp1; #Replicate the simulated coefficient J times, SxJ
    β_1_t_temp = β_1_t.*temp1; #Replicate the simulated coefficient J times, SxJ
    β_2_t_temp = β_2_t.*temp1; #Replicate the simulated coefficient J times, SxJ
    β_3_t_temp = abs.(β_3_t.*temp1); #Replicate the simulated coefficient J times, SxJ
    ξ_t_temp =  charProduct_2[1,:,t]'.*temp1; # SxJ


    inContrNum = exp.(β_0_t_temp+ β_1_t_temp.* charProduct_1[1,:,t]'+β_2_t_temp.* charProduct_1[2,:,t]'-β_3_t_temp.*price'+ξ_t_temp);
    inContrNum = hcat(ones(S), inContrNum)

    total_sum = sum.(inContrNum[s,:] for s in 1:S);
    f_jt = inContrNum./total_sum; #SxJ
    sum_f_jt = sum.(f_jt[:,p] for p in 1:J+1)
    mShare = (1/S)*sum_f_jt

    derivative = f_jt .* (ones(S,J+1) - f_jt) .* (-1) .* (ones(S,J+1) .* β_3_t) 
    sum_derivative = sum.(derivative[:,p] for p in 1:J+1)

    finDerivative   =   (1/S) .* sum_derivative   
    
    return mShare, finDerivative
end