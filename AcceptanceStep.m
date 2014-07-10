function [BoolAccept] = AcceptanceStep(iter,Fsimvalues,ETA,Evaluated_Points,LastAcceptedPoint,Betas)

temp = Evaluated_Points(1,:);
PROBLEMDIMENSION = length(temp);

function [z] = obj(x0)
    z= Betas(1) + Betas(2 : PROBLEMDIMENSION+1)'*x0 + Betas(PROBLEMDIMENSION + 2 : 2*PROBLEMDIMENSION+1)'*x0.^2;

end

if(iter == 1)
    BoolAccept=1;
else
    
     
    PredictedReduction = obj(Evaluated_Points(LastAcceptedPoint,:)') - obj(Evaluated_Points(iter,:)');
    ActualReduction = Fsimvalues(LastAcceptedPoint) - Fsimvalues(iter);
    
    Ratio = ActualReduction/PredictedReduction;
    if (Ratio > ETA && ActualReduction > 0)
        BoolAccept=1;
    else
        BoolAccept=0;
    end
    
end

end
