function wucha=calculation_error(jieguo,test,h); 
global wuchaid;
%%% id=[id expfloor];
[mx,my]=size(jieguo);%%%%%minx,miny,maxx,maxy
[hang,lie]=size(test);
if (mx==hang)&(my==lie)
	if wuchaid==1
            wucha=sqrt(sum(sum((jieguo-test).^2))/(mx*my));
        elseif wuchaid==2
            wucha=(sum(sum(abs(jieguo-test)))/(mx*my));  
        else wuchaid==3
            wucha=(sum(sum((jieguo-test)))/(mx*my));  
	end
else
    temp=0;
    x2=round2(test(:,1)/h)+1;%%%%%%%%%
    y2=round2(test(:,2)/h)+1;
    for i=1:hang
            temp=jieguo(x2(i),y2(i));
            wuchas(i)=(temp-test(i,3));
    end
    if wuchaid==1
        wucha=sqrt(sum(wuchas.^2)/hang);
    elseif wuchaid==2
        wucha=sum(abs(wuchas))/hang;
    else wuchaid==3
        wucha=sum(wuchas)/hang;
    end
end