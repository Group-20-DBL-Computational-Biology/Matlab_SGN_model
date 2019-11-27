 function dydt = SystemState(t,y, para)
       
        
        dydt = [
         0;
        ((para(1) * y(5)) * y(1))/(para(7) + y(1)) - (para(4) * y(2));
        (para(2) * y(6) * y(2))/(para(8) + y(2))- (para(6) * y(3));
        para(6) * y(3);
        -(para(3) * y(5) * y(1))/(para(7) + y(1));
        -(para(5) * y(6)) /(para(9) + y(6))] ;        
        
       
end 