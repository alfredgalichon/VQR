function [] = PlotFixedx2D( xeval,beta,lT,save,namefile,pref, suff,suff2) 
yeval=squeeze(yhat2D(xeval,beta));
titlesize=16
figure
surf(reshape(yeval(:,1),lT-1,lT-1));
xlabel('u_1')
ylabel('u_2')
zlabel('y_1')
name= '$(u_1, u_2) \mapsto \beta_1(u_1, u_2)^\top x$';
%name=strcat('$(u_1,u_2) \to ',pref,'_{1',suff,suff2);
title(name, 'Interpreter','latex','fontsize',titlesize)
if save
    h=gcf;
    print(h,'-dpdf',strcat(namefile,'1',suff,'.pdf'))
end
figure
surf(reshape(yeval(:,2),lT-1,lT-1));
xlabel('u_1')
ylabel('u_2')
zlabel('y_2')
%name=strcat('$(u_1,u_2) \to ',pref,'_{2',suff,suff2);
name= '$(u_1, u_2) \mapsto \beta_2(u_1, u_2)^\top x$';
title(name, 'Interpreter','latex','fontsize',titlesize)
if save
    h=gcf;
    print(h,'-dpdf',strcat(namefile,'2',suff,'.pdf'))
end

end

