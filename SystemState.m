 function dydt = SystemState(t,y, para)
        %% parameters
        %link each concentration to the right parameter
        DNA = y(1);
        mRNA_lacI = y(2);
        mRNA_lacZ = y(3);
        TsR = y(4);
        TlR = y(5);
        RA = y(6);
        R = y(7); 
        O = y(8);
        RO = y(9);
        L = y(10);
        G = y(11); 
        GL = y(12); 
        A = y(13);
        
        %link each parameter value 
        kts = para(1); %nanomolar per minute
        ktl = para(2); %nanomolar per minute
        kcs = para(3); % per minute
        deltamRNA = para(4); %per minute 
        deltaTlr = para(5); %per minute
        Ks = para(6);% nanomolar
        Kl = para(7); %nanomolar
        Ktlr = para(8); %
        
        %degradation values
        deltaR = para(9); %not known yet
        deltaG = para(10); 
        
        %k-values for repressor and allolactose
        k_on1 = para(11);
        k_off1 = para(12);
        
        %k-values repressor and operator
        k_on2 = para(13);
        k_off2 = para(14);
        
        %k-values galactosidase and lactose
        k_on3 = para(15); 
        k_off3 = para(16); 
        kcat3 = para(17); 
        %% ODEs
        dDNA = 0; 
        %ODEs needed for the expression lacI
        %transcription
        dmRNA_lacI = (kts * TsR * DNA)/(Ks + DNA) - (deltamRNA * mRNA_lacI);
        %translation + degradation
        trans_R = (ktl * TlR * mRNA_lacI)/(Kl + mRNA_lacI) -deltaR* R;
              
        %ODEs needed for the expression lacZ
        %transcription
        dmRNA_lacZ = (kts * TsR * O)/(Ks + O) - (deltamRNA * mRNA_lacZ);
        %translation  + degradation
        trans_G = (ktl * TlR * mRNA_lacZ)/(Kl + mRNA_lacZ)- deltaG * G;
        
        %Resources for both expressions
        dTsR = -(kcs * TsR * DNA)/(Ks + DNA);
        dTlR = -(deltaTlr * TlR) /(Ktlr + TlR);
        
        %1) ODEs neede for itneraction between repressor and allolactose 
        rep_no_allo = -k_on1 * R * A + k_off1 * RA; 
        free_allo = -k_on1 * R * A + k_off1 * RA; 
        dRA = k_on1 * R * A -k_off1 * RA; 
        
        %2) ODEs needed for interaction between repressor and operator
        dR = (-k_on2 * R * O + k_off2 * RO) + trans_R + rep_no_allo; %free repressor
        dO = (-k_on2 * R * O + k_off2 * RO); %free DNA = free operator
        dRO = (k_on2 * R * O + k_off2 *RO); %complex
        
        %ODEs needed for interactions between galactosidase and lactose
        dL = -k_on3*G * L + k_off3* GL;
        dG = trans_G + (-k_on3 * G * L + (k_off3 + kcat3) * GL);
        dGL = k_on3*G * L - (k_off3 + kcat3) * GL; 
        dA = 0.5*kcat3 * GL + free_allo; 

        %store all solutions in the following array 
        dydt = [dDNA; dmRNA_lacI; dmRNA_lacZ; dTsR; dTlR; dRA; dR; dO; dRO; dL; dG; dGL; dA];
end 