function lisangzhi=climate_site_discrete degrees(x,y,z);%%%%Calculate the discrete values of the points in the neighborhood of each point based on the delaunay Tyson triangulation
%%%% lamdaid=0 Indicates that all sampling points have the same weight
global lamdazhi lamdaid;
R=min(15,length(x)/2);
m=length(x);
if lamdaid==0
    lisangzhi=lamdazhi*ones(m,1);
elseif lamdaid==1%%%Determine the weight of this point based on the variance of this sample point and its linked sample points 
    tri=delaunay(x,y);
    for i=1:m
        tt1=find(tri(:,1)==i);
        tt2=find(tri(:,2)==i);
        tt3=find(tri(:,3)==i);
        dd=union(union(tt1,tt2),tt3);  %%%Row numbers of all triangles tri containing i
        ddi=union(union(tri(dd,1),tri(dd,2)),tri(dd,3));   %%%% All elements of tri rows containing i
        ddi=setdiff(ddi,i);  %%%Remove duplicates from
        lisangzhi(i)=var(z(ddi)); %%%%variance (statistics)
        clear dd ddi;
    end
    lisangzhi=lisangzhi'; %%%%%  The output is a column vector
    maxzhi=max(lisangzhi);
    minzhi=min(lisangzhi);
    if maxzhi==minzhi
        lisangzhi=lamdazhi*ones(length(lisangzhi),1);
    else
        %lisangzhi=(lisangzhi*(99/(maxzhi-minzhi))+(maxzhi-100*minzhi)/(maxzhi-minzhi))*lamda/100;
        lisangzhi=(-lisangzhi*(99/(maxzhi-minzhi))+(100*maxzhi-minzhi)/(maxzhi-minzhi))*lamdazhi/100;%%%%The largest variance is 0.01 and the smallest variance is 1
    end
elseif lamdaid==2%%The weight of this point is determined based on the difference between the mean values of this sample point and its surrounding neighbors, with the aim of removing the effect of certain anomalies
    for i=1:m
        juli=(x-x(i)).^2+(y-y(i)).^2;
        [sortjuli,ID]=sort(juli);    %%%Distance Sort
        newz=z(ID(2:R));
        cha(i)=abs(z(i)-mean(newz));
        clear newz juli sortjuli ID;
    end
    maxzhi=max(cha);
    minzhi=min(cha);
    if maxzhi==minzhi
        lisangzhi=lamdazhi*ones(length(x),1);
    else
        lisangzhi=lamdazhi*(maxzhi-0.00001*minzhi-cha*0.99999)/(maxzhi-minzhi);%%%%The largest difference, 0.00001, and the smallest difference, 1
    end
end