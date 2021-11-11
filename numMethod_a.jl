function FOC1(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10)
    s  = 0.0 #Initialize
    ds = 0.0 #Initialize
    β_0_t = β_0[:,t];
    β_1_t = β_1[:,t];
    β_2_t = β_2[:,t];
    β_3_t = β_3[:,t];
    ξ_t = charProduct_2[1,:,t];

    for i in 1:S

        prod1 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,1,t] + β_2_t[i]*charProduct_1[2,1,t] - (abs(β_3_t[i])*p_1) + ξ_t[1])
        prod2 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,2,t] + β_2_t[i]*charProduct_1[2,2,t] - (abs(β_3_t[i])*p_2) + ξ_t[2])
        prod3 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,3,t] + β_2_t[i]*charProduct_1[2,3,t] - (abs(β_3_t[i])*p_3) + ξ_t[3])
        prod4 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,4,t] + β_2_t[i]*charProduct_1[2,4,t] - (abs(β_3_t[i])*p_4) + ξ_t[4])
        prod5 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,5,t] + β_2_t[i]*charProduct_1[2,5,t] - (abs(β_3_t[i])*p_5) + ξ_t[5])
        prod6 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,6,t] + β_2_t[i]*charProduct_1[2,6,t] - (abs(β_3_t[i])*p_6) + ξ_t[6])
        prod7 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,7,t] + β_2_t[i]*charProduct_1[2,7,t] - (abs(β_3_t[i])*p_7) + ξ_t[7])
        prod8 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,8,t] + β_2_t[i]*charProduct_1[2,8,t] - (abs(β_3_t[i])*p_8) + ξ_t[8])
        prod9 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,9,t] + β_2_t[i]*charProduct_1[2,9,t] - (abs(β_3_t[i])*p_9) + ξ_t[9])
        prod10 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,10,t] + β_2_t[i]*charProduct_1[2,10,t] - (abs(β_3_t[i])*p_10) + ξ_t[10])

        denom = 1 + prod1 + prod2 + prod3 + prod4 + prod5 + prod6 + prod7 + prod8 + prod9 + prod10 #denominator for f_jt

        ds = ds - (1/S)*(prod1/denom)*(1-prod1/denom)*abs(β_3_t[i]) #Derivative of s with respect to p
        s  = s + (1/S)*prod1/denom; # Sum of f_jt with respect to individual i

    end

    output1  = s + (p_1 - mc[1,t])*ds

    return output1

end


function FOC2(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10)
    s  = 0.0 #Initialize
    ds = 0.0 #Initialize

    β_0_t = β_0[:,t];
    β_1_t = β_1[:,t];
    β_2_t = β_2[:,t];
    β_3_t = β_3[:,t];
    ξ_t = charProduct_2[1,:,t];

    for i in 1:S

        prod1 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,1,t] + β_2_t[i]*charProduct_1[2,1,t] - (abs(β_3_t[i])*p_1) + ξ_t[1])
        prod2 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,2,t] + β_2_t[i]*charProduct_1[2,2,t] - (abs(β_3_t[i])*p_2) + ξ_t[2])
        prod3 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,3,t] + β_2_t[i]*charProduct_1[2,3,t] - (abs(β_3_t[i])*p_3) + ξ_t[3])
        prod4 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,4,t] + β_2_t[i]*charProduct_1[2,4,t] - (abs(β_3_t[i])*p_4) + ξ_t[4])
        prod5 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,5,t] + β_2_t[i]*charProduct_1[2,5,t] - (abs(β_3_t[i])*p_5) + ξ_t[5])
        prod6 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,6,t] + β_2_t[i]*charProduct_1[2,6,t] - (abs(β_3_t[i])*p_6) + ξ_t[6])
        prod7 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,7,t] + β_2_t[i]*charProduct_1[2,7,t] - (abs(β_3_t[i])*p_7) + ξ_t[7])
        prod8 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,8,t] + β_2_t[i]*charProduct_1[2,8,t] - (abs(β_3_t[i])*p_8) + ξ_t[8])
        prod9 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,9,t] + β_2_t[i]*charProduct_1[2,9,t] - (abs(β_3_t[i])*p_9) + ξ_t[9])
        prod10 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,10,t] + β_2_t[i]*charProduct_1[2,10,t] - (abs(β_3_t[i])*p_10) + ξ_t[10])

        denom = 1 + prod1 + prod2 + prod3 + prod4 + prod5 + prod6 + prod7 + prod8 + prod9 + prod10 #denominator for f_jt

        ds = ds - (1/S)*(prod2/denom)*(1-prod2/denom)*abs(β_3_t[i]) #Derivative of s with respect to p
        s  = s + (1/S)*prod2/denom; # Sum of f_jt with respect to individual i

    end

    output2  = s + (p_2 - mc[2,t])*ds

    return output2

end

function FOC3(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10)
    s  = 0.0 #Initialize
    ds = 0.0 #Initialize

    β_0_t = β_0[:,t];
    β_1_t = β_1[:,t];
    β_2_t = β_2[:,t];
    β_3_t = β_3[:,t];
    ξ_t = charProduct_2[1,:,t];

    for i in 1:S

        prod1 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,1,t] + β_2_t[i]*charProduct_1[2,1,t] - (abs(β_3_t[i])*p_1) + ξ_t[1])
        prod2 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,2,t] + β_2_t[i]*charProduct_1[2,2,t] - (abs(β_3_t[i])*p_2) + ξ_t[2])
        prod3 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,3,t] + β_2_t[i]*charProduct_1[2,3,t] - (abs(β_3_t[i])*p_3) + ξ_t[3])
        prod4 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,4,t] + β_2_t[i]*charProduct_1[2,4,t] - (abs(β_3_t[i])*p_4) + ξ_t[4])
        prod5 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,5,t] + β_2_t[i]*charProduct_1[2,5,t] - (abs(β_3_t[i])*p_5) + ξ_t[5])
        prod6 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,6,t] + β_2_t[i]*charProduct_1[2,6,t] - (abs(β_3_t[i])*p_6) + ξ_t[6])
        prod7 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,7,t] + β_2_t[i]*charProduct_1[2,7,t] - (abs(β_3_t[i])*p_7) + ξ_t[7])
        prod8 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,8,t] + β_2_t[i]*charProduct_1[2,8,t] - (abs(β_3_t[i])*p_8) + ξ_t[8])
        prod9 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,9,t] + β_2_t[i]*charProduct_1[2,9,t] - (abs(β_3_t[i])*p_9) + ξ_t[9])
        prod10 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,10,t] + β_2_t[i]*charProduct_1[2,10,t] - (abs(β_3_t[i])*p_10) + ξ_t[10])

        denom = 1 + prod1 + prod2 + prod3 + prod4 + prod5 + prod6 + prod7 + prod8 + prod9 + prod10 #denominator for f_jt

        ds = ds - (1/S)*(prod3/denom)*(1-prod3/denom)*abs(β_3_t[i]) #Derivative of s with respect to p
        s  = s + (1/S)*prod3/denom; # Sum of f_jt with respect to individual i

    end

    output3  = s + (p_3 - mc[3,t])*ds

    return output3

end

function FOC4(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10)
    s  = 0.0 #Initialize
    ds = 0.0 #Initialize

    β_0_t = β_0[:,t];
    β_1_t = β_1[:,t];
    β_2_t = β_2[:,t];
    β_3_t = β_3[:,t];
    ξ_t = charProduct_2[1,:,t];

    for i in 1:S

        prod1 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,1,t] + β_2_t[i]*charProduct_1[2,1,t] - (abs(β_3_t[i])*p_1) + ξ_t[1])
        prod2 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,2,t] + β_2_t[i]*charProduct_1[2,2,t] - (abs(β_3_t[i])*p_2) + ξ_t[2])
        prod3 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,3,t] + β_2_t[i]*charProduct_1[2,3,t] - (abs(β_3_t[i])*p_3) + ξ_t[3])
        prod4 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,4,t] + β_2_t[i]*charProduct_1[2,4,t] - (abs(β_3_t[i])*p_4) + ξ_t[4])
        prod5 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,5,t] + β_2_t[i]*charProduct_1[2,5,t] - (abs(β_3_t[i])*p_5) + ξ_t[5])
        prod6 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,6,t] + β_2_t[i]*charProduct_1[2,6,t] - (abs(β_3_t[i])*p_6) + ξ_t[6])
        prod7 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,7,t] + β_2_t[i]*charProduct_1[2,7,t] - (abs(β_3_t[i])*p_7) + ξ_t[7])
        prod8 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,8,t] + β_2_t[i]*charProduct_1[2,8,t] - (abs(β_3_t[i])*p_8) + ξ_t[8])
        prod9 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,9,t] + β_2_t[i]*charProduct_1[2,9,t] - (abs(β_3_t[i])*p_9) + ξ_t[9])
        prod10 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,10,t] + β_2_t[i]*charProduct_1[2,10,t] - (abs(β_3_t[i])*p_10) + ξ_t[10])

        denom = 1 + prod1 + prod2 + prod3 + prod4 + prod5 + prod6 + prod7 + prod8 + prod9 + prod10 #denominator for f_jt

        ds = ds - (1/S)*(prod4/denom)*(1-prod4/denom)*abs(β_3_t[i]) #Derivative of s with respect to p
        s  = s + (1/S)*prod4/denom; # Sum of f_jt with respect to individual i

    end

    output4  = s + (p_4 - mc[4,t])*ds

    return output4

end

function FOC5(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10)
    s  = 0.0 #Initialize
    ds = 0.0 #Initialize

    β_0_t = β_0[:,t];
    β_1_t = β_1[:,t];
    β_2_t = β_2[:,t];
    β_3_t = β_3[:,t];
    ξ_t = charProduct_2[1,:,t];

    for i in 1:S

        prod1 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,1,t] + β_2_t[i]*charProduct_1[2,1,t] - (abs(β_3_t[i])*p_1) + ξ_t[1])
        prod2 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,2,t] + β_2_t[i]*charProduct_1[2,2,t] - (abs(β_3_t[i])*p_2) + ξ_t[2])
        prod3 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,3,t] + β_2_t[i]*charProduct_1[2,3,t] - (abs(β_3_t[i])*p_3) + ξ_t[3])
        prod4 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,4,t] + β_2_t[i]*charProduct_1[2,4,t] - (abs(β_3_t[i])*p_4) + ξ_t[4])
        prod5 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,5,t] + β_2_t[i]*charProduct_1[2,5,t] - (abs(β_3_t[i])*p_5) + ξ_t[5])
        prod6 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,6,t] + β_2_t[i]*charProduct_1[2,6,t] - (abs(β_3_t[i])*p_6) + ξ_t[6])
        prod7 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,7,t] + β_2_t[i]*charProduct_1[2,7,t] - (abs(β_3_t[i])*p_7) + ξ_t[7])
        prod8 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,8,t] + β_2_t[i]*charProduct_1[2,8,t] - (abs(β_3_t[i])*p_8) + ξ_t[8])
        prod9 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,9,t] + β_2_t[i]*charProduct_1[2,9,t] - (abs(β_3_t[i])*p_9) + ξ_t[9])
        prod10 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,10,t] + β_2_t[i]*charProduct_1[2,10,t] - (abs(β_3_t[i])*p_10) + ξ_t[10])

        denom = 1 + prod1 + prod2 + prod3 + prod4 + prod5 + prod6 + prod7 + prod8 + prod9 + prod10 #denominator for f_jt

        ds = ds - (1/S)*(prod5/denom)*(1-prod5/denom)*abs(β_3_t[i]) #Derivative of s with respect to p
        s  = s + (1/S)*prod5/denom; # Sum of f_jt with respect to individual i

    end

    output5  = s + (p_5 - mc[5,t])*ds

    return output5

end

function FOC6(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10)
    s  = 0.0 #Initialize
    ds = 0.0 #Initialize

    β_0_t = β_0[:,t];
    β_1_t = β_1[:,t];
    β_2_t = β_2[:,t];
    β_3_t = β_3[:,t];
    ξ_t = charProduct_2[1,:,t];

   for i in 1:S

        prod1 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,1,t] + β_2_t[i]*charProduct_1[2,1,t] - (abs(β_3_t[i])*p_1) + ξ_t[1])
        prod2 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,2,t] + β_2_t[i]*charProduct_1[2,2,t] - (abs(β_3_t[i])*p_2) + ξ_t[2])
        prod3 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,3,t] + β_2_t[i]*charProduct_1[2,3,t] - (abs(β_3_t[i])*p_3) + ξ_t[3])
        prod4 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,4,t] + β_2_t[i]*charProduct_1[2,4,t] - (abs(β_3_t[i])*p_4) + ξ_t[4])
        prod5 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,5,t] + β_2_t[i]*charProduct_1[2,5,t] - (abs(β_3_t[i])*p_5) + ξ_t[5])
        prod6 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,6,t] + β_2_t[i]*charProduct_1[2,6,t] - (abs(β_3_t[i])*p_6) + ξ_t[6])
        prod7 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,7,t] + β_2_t[i]*charProduct_1[2,7,t] - (abs(β_3_t[i])*p_7) + ξ_t[7])
        prod8 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,8,t] + β_2_t[i]*charProduct_1[2,8,t] - (abs(β_3_t[i])*p_8) + ξ_t[8])
        prod9 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,9,t] + β_2_t[i]*charProduct_1[2,9,t] - (abs(β_3_t[i])*p_9) + ξ_t[9])
        prod10 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,10,t] + β_2_t[i]*charProduct_1[2,10,t] - (abs(β_3_t[i])*p_10) + ξ_t[10])

        denom = 1 + prod1 + prod2 + prod3 + prod4 + prod5 + prod6 + prod7 + prod8 + prod9 + prod10 #denominator for f_jt

        ds = ds - (1/S)*(prod6/denom)*(1-prod6/denom)*abs(β_3_t[i]) #Derivative of s with respect to p
        s  = s + (1/S)*prod6/denom; # Sum of f_jt with respect to individual i

    end

    output6  = s + (p_6 - mc[6,t])*ds

    return output6

end

function FOC7(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10)
    s  = 0.0 #Initialize
    ds = 0.0 #Initialize

    β_0_t = β_0[:,t];
    β_1_t = β_1[:,t];
    β_2_t = β_2[:,t];
    β_3_t = β_3[:,t];
    ξ_t = charProduct_2[1,:,t];

    for i in 1:S

        prod1 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,1,t] + β_2_t[i]*charProduct_1[2,1,t] - (abs(β_3_t[i])*p_1) + ξ_t[1])
        prod2 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,2,t] + β_2_t[i]*charProduct_1[2,2,t] - (abs(β_3_t[i])*p_2) + ξ_t[2])
        prod3 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,3,t] + β_2_t[i]*charProduct_1[2,3,t] - (abs(β_3_t[i])*p_3) + ξ_t[3])
        prod4 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,4,t] + β_2_t[i]*charProduct_1[2,4,t] - (abs(β_3_t[i])*p_4) + ξ_t[4])
        prod5 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,5,t] + β_2_t[i]*charProduct_1[2,5,t] - (abs(β_3_t[i])*p_5) + ξ_t[5])
        prod6 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,6,t] + β_2_t[i]*charProduct_1[2,6,t] - (abs(β_3_t[i])*p_6) + ξ_t[6])
        prod7 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,7,t] + β_2_t[i]*charProduct_1[2,7,t] - (abs(β_3_t[i])*p_7) + ξ_t[7])
        prod8 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,8,t] + β_2_t[i]*charProduct_1[2,8,t] - (abs(β_3_t[i])*p_8) + ξ_t[8])
        prod9 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,9,t] + β_2_t[i]*charProduct_1[2,9,t] - (abs(β_3_t[i])*p_9) + ξ_t[9])
        prod10 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,10,t] + β_2_t[i]*charProduct_1[2,10,t] - (abs(β_3_t[i])*p_10) + ξ_t[10])

        denom = 1 + prod1 + prod2 + prod3 + prod4 + prod5 + prod6 + prod7 + prod8 + prod9 + prod10 #denominator for f_jt

        ds = ds - (1/S)*(prod7/denom)*(1-prod7/denom)*abs(β_3_t[i]) #Derivative of s with respect to p
        s  = s + (1/S)*prod7/denom; # Sum of f_jt with respect to individual i

    end

    output7  = s + (p_7 - mc[7,t])*ds

    return output7

end

function FOC8(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10)
    s  = 0.0 #Initialize
    ds = 0.0 #Initialize

    β_0_t = β_0[:,t];
    β_1_t = β_1[:,t];
    β_2_t = β_2[:,t];
    β_3_t = β_3[:,t];
    ξ_t = charProduct_2[1,:,t];

    for i in 1:S

        prod1 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,1,t] + β_2_t[i]*charProduct_1[2,1,t] - (abs(β_3_t[i])*p_1) + ξ_t[1])
        prod2 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,2,t] + β_2_t[i]*charProduct_1[2,2,t] - (abs(β_3_t[i])*p_2) + ξ_t[2])
        prod3 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,3,t] + β_2_t[i]*charProduct_1[2,3,t] - (abs(β_3_t[i])*p_3) + ξ_t[3])
        prod4 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,4,t] + β_2_t[i]*charProduct_1[2,4,t] - (abs(β_3_t[i])*p_4) + ξ_t[4])
        prod5 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,5,t] + β_2_t[i]*charProduct_1[2,5,t] - (abs(β_3_t[i])*p_5) + ξ_t[5])
        prod6 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,6,t] + β_2_t[i]*charProduct_1[2,6,t] - (abs(β_3_t[i])*p_6) + ξ_t[6])
        prod7 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,7,t] + β_2_t[i]*charProduct_1[2,7,t] - (abs(β_3_t[i])*p_7) + ξ_t[7])
        prod8 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,8,t] + β_2_t[i]*charProduct_1[2,8,t] - (abs(β_3_t[i])*p_8) + ξ_t[8])
        prod9 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,9,t] + β_2_t[i]*charProduct_1[2,9,t] - (abs(β_3_t[i])*p_9) + ξ_t[9])
        prod10 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,10,t] + β_2_t[i]*charProduct_1[2,10,t] - (abs(β_3_t[i])*p_10) + ξ_t[10])

        denom = 1 + prod1 + prod2 + prod3 + prod4 + prod5 + prod6 + prod7 + prod8 + prod9 + prod10 #denominator for f_jt

        ds = ds - (1/S)*(prod8/denom)*(1-prod8/denom)*abs(β_3_t[i]) #Derivative of s with respect to p
        s  = s + (1/S)*prod8/denom; # Sum of f_jt with respect to individual i

    end

    output8  = s + (p_8 - mc[8,t])*ds

    return output8

end

function FOC9(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10)
    s  = 0.0 #Initialize
    ds = 0.0 #Initialize

    β_0_t = β_0[:,t];
    β_1_t = β_1[:,t];
    β_2_t = β_2[:,t];
    β_3_t = β_3[:,t];
    ξ_t = charProduct_2[1,:,t];

    for i in 1:S

        prod1 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,1,t] + β_2_t[i]*charProduct_1[2,1,t] - (abs(β_3_t[i])*p_1) + ξ_t[1])
        prod2 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,2,t] + β_2_t[i]*charProduct_1[2,2,t] - (abs(β_3_t[i])*p_2) + ξ_t[2])
        prod3 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,3,t] + β_2_t[i]*charProduct_1[2,3,t] - (abs(β_3_t[i])*p_3) + ξ_t[3])
        prod4 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,4,t] + β_2_t[i]*charProduct_1[2,4,t] - (abs(β_3_t[i])*p_4) + ξ_t[4])
        prod5 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,5,t] + β_2_t[i]*charProduct_1[2,5,t] - (abs(β_3_t[i])*p_5) + ξ_t[5])
        prod6 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,6,t] + β_2_t[i]*charProduct_1[2,6,t] - (abs(β_3_t[i])*p_6) + ξ_t[6])
        prod7 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,7,t] + β_2_t[i]*charProduct_1[2,7,t] - (abs(β_3_t[i])*p_7) + ξ_t[7])
        prod8 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,8,t] + β_2_t[i]*charProduct_1[2,8,t] - (abs(β_3_t[i])*p_8) + ξ_t[8])
        prod9 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,9,t] + β_2_t[i]*charProduct_1[2,9,t] - (abs(β_3_t[i])*p_9) + ξ_t[9])
        prod10 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,10,t] + β_2_t[i]*charProduct_1[2,10,t] - (abs(β_3_t[i])*p_10) + ξ_t[10])

        denom = 1 + prod1 + prod2 + prod3 + prod4 + prod5 + prod6 + prod7 + prod8 + prod9 + prod10 #denominator for f_jt

        ds = ds - (1/S)*(prod9/denom)*(1-prod9/denom)*abs(β_3_t[i]) #Derivative of s with respect to p
        s  = s + (1/S)*prod9/denom; # Sum of f_jt with respect to individual i

    end

    output9  = s + (p_9 - mc[9,t])*ds

    return output9

end

function FOC10(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10)
    s  = 0.0 #Initialize
    ds = 0.0 #Initialize

    β_0_t = β_0[:,t];
    β_1_t = β_1[:,t];
    β_2_t = β_2[:,t];
    β_3_t = β_3[:,t];
    ξ_t = charProduct_2[1,:,t];

    for i in 1:S

        prod1 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,1,t] + β_2_t[i]*charProduct_1[2,1,t] - (abs(β_3_t[i])*p_1) + ξ_t[1])
        prod2 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,2,t] + β_2_t[i]*charProduct_1[2,2,t] - (abs(β_3_t[i])*p_2) + ξ_t[2])
        prod3 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,3,t] + β_2_t[i]*charProduct_1[2,3,t] - (abs(β_3_t[i])*p_3) + ξ_t[3])
        prod4 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,4,t] + β_2_t[i]*charProduct_1[2,4,t] - (abs(β_3_t[i])*p_4) + ξ_t[4])
        prod5 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,5,t] + β_2_t[i]*charProduct_1[2,5,t] - (abs(β_3_t[i])*p_5) + ξ_t[5])
        prod6 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,6,t] + β_2_t[i]*charProduct_1[2,6,t] - (abs(β_3_t[i])*p_6) + ξ_t[6])
        prod7 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,7,t] + β_2_t[i]*charProduct_1[2,7,t] - (abs(β_3_t[i])*p_7) + ξ_t[7])
        prod8 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,8,t] + β_2_t[i]*charProduct_1[2,8,t] - (abs(β_3_t[i])*p_8) + ξ_t[8])
        prod9 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,9,t] + β_2_t[i]*charProduct_1[2,9,t] - (abs(β_3_t[i])*p_9) + ξ_t[9])
        prod10 = exp(β_0_t[i] + β_1_t[i]*charProduct_1[1,10,t] + β_2_t[i]*charProduct_1[2,10,t] - (abs(β_3_t[i])*p_10) + ξ_t[10])

        denom = 1 + prod1 + prod2 + prod3 + prod4 + prod5 + prod6 + prod7 + prod8 + prod9 + prod10 #denominator for f_jt

        ds = ds - (1/S)*(prod10/denom)*(1-prod10/denom)*abs(β_3_t[i]) #Derivative of s with respect to p
        s  = s + (1/S)*prod10/denom; # Sum of f_jt with respect to individual i

    end

    output10  = s + (p_10 - mc[10,t])*ds

    return output10

end