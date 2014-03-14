function [BoolAccept] = AcceptanceStep(iter,Fsimvalues,ETA,Evaluated_Points,LastAcceptedPoint,Betas)

function [z] = obj(x0)
    z=0;
    z = z+ Betas'*x0;
end

if(iter == 1)
    BoolAccept=1;
else
    
     
    PredictedReduction = obj(Evaluated_Points(LastAcceptedPoint,:)') - obj(Evaluated_Points(iter,:)');
    ActualReduction = Fsimvalues(LastAcceptedPoint) - Fsimvalues(iter);
    
    Ratio = ActualReduction/PredictedReduction;
    if (Ratio > ETA)
        BoolAccept=1;
    else
        BoolAccept=0;
    end
    
end

end
