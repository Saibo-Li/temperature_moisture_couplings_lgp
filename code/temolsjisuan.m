clear all;  %%%%利用IDW和KRI残差计算最小二乘方法+残差插值获得的江西地表气温
mulu='C:\duzp\guokeda\jxtemdata\';
demmulu=strcat(mulu,'jxexdem.txt');
krimulu=strcat(mulu,'cckri.txt');
idwmulu=strcat(mulu,'ccidw.txt');
[X0,Y0,H,mx,my,dem]=shuju2(demmulu);
[X0,Y0,H,mx,my,cckri]=shuju2(krimulu);
[X0,Y0,H,mx,my,ccidw]=shuju2(idwmulu);
idwnew=-9999*ones(mx,my);
krinew=-9999*ones(mx,my);
for i=1:mx
    for j=1:my
        xx=X0+(i-1)*H+0.5*H;
        yy=Y0+(j-1)*H+0.5*H;
        if dem(i,j)>-9999
            idwnew(i,j)=19.590972+xx*0.000002-yy*0.000007-dem(i,j)*0.004563+ccidw(i,j);
            krinew(i,j)=19.590972+xx*0.000002-yy*0.000007-dem(i,j)*0.004563+cckri(i,j);
        end
    end
end

jieguoshuchu(idwnew,strcat(mulu,'jxtemidw.txt'),mx,my,X0,Y0,H);
jieguoshuchu(krinew,strcat(mulu,'jxtemkri.txt'),mx,my,X0,Y0,H);
