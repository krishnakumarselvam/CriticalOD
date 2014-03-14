function BoolMixturePoint = EvaluateChangeinBeta(OldBeta,CurrBeta,BETATOL)

    relBeta = norm(CurrBeta-OldBeta)/norm(OldBeta);
    if (relBeta < BETATOL)
        BoolMixturePoint =1;
    else
        BoolMixturePoint=0;
    end
end