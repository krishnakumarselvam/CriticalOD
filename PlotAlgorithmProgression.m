function [] = PlotAlgorithmProgression(Fsimvalues_Secondary,Fsimvalues,AcceptedInSO,iter,AcceptedFsimValues,IsMixturePoint)
  close all
 
  
  grey = [0.4 , 0.4 , 0.4] ;
  plot(Fsimvalues(1:iter),'-*','color',grey);
  hold on
  plot(AcceptedFsimValues(1:iter),'-*','color','blue');
  loc = find(IsMixturePoint);
  plot(loc,Fsimvalues(loc),'r*');
 
  title([num2str(sum(AcceptedInSO)) ' points accepted out of ' num2str(iter) ' samples']);
  xlabel('Iterations');
  ylabel('Simulated Values');
  legend('Simulated Values (SDev)', 'Accepted Values (SDev)', 'Mixture Points (SDev)');
  
  figure 
  plot(Fsimvalues_Secondary(1:iter),'-*','color','black');
  xlabel('Iterations');
  ylabel('Simulated Values');
  legend ('Average Travel Time');
    
end