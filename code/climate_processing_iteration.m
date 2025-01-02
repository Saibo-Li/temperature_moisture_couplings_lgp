function jieguo=hasmguding(chazhi,caiyang,h);
global caiyangid guimo hasmtime yuchuliid luid quanzhong;
%%%%  HASM iterative program with fixed number of iterations and no error test
[mx,my]=size(chazhi);
if isempty(caiyang(:,1))
    disp('There are no sampling points in the area, check the program and data for errors!');
    jieguo=chazhi;
else
    lamda=lisangdu(caiyang(:,1),caiyang(:,2),caiyang(:,3));%%%%Compute lambda vectors
    if ismember(caiyangid,[1 3 5])==0
        disp('caiyangid assignment error (can only be 1, 3 or 5), the program automatically changed to 3');
        caiyangid=3;
    end
    [GG,GH]=taylorcaiyang(caiyang(:,1),caiyang(:,2),caiyang(:,3),chazhi,lamda,h);%%%%Compute the sampling vector and sampling matrix to generate lambda^2*G'G,lambda^2*G'*H
   
    if (mx*my)<guimo  %%%Small scale, use direct solution
        xishus=taylorxishu(mx,my,GG);%%%%Coefficient matrix A'*A+B'*B+lambda^2*G'G
        clear GG;
%         [l,u]=shangxiajie3(caiyang,grid,mx,my,h);
        if luid==1  %%%Are upper and lower boundary controls used?
            [l,u]=shangxiajie(caiyang,mx,my,h);%%%%Generate upper and lower bounds
        end
         clear caiyang grid;
        for is=1:hasmtime
            b=quyumoni(chazhi,h)+GH;  %%%%Right end term A'*c + B'*d + lambda^2*G'*H
            xx=xishus\b;
            clear b;
            chazhi(2:(mx-1),2:(my-1))=(reshape(xx,(my-2),(mx-2)))';%%%%Take a one-dimensional vector and turn it into a two-dimensional matrix.
            clear xx;
            if luid==1
                for i=1:(mx-2)
                    for j=1:(my-2)
                        if chazhi(i+1,j+1)<l(i,j)
                            chazhi(i+1,j+1)=l(i,j);
                        elseif chazhi(i+1,j+1)>u(i,j)
                            chazhi(i+1,j+1)=u(i,j);
                        end
                    end
                end
            end
        end
    else %%%Large scale with PCG solution
        M=duijiaoxian(mx,my); %%% Main diagonal of A'A + B'B
        yuxishu=yuchulixishu(mx,my,GG,M,yuchuliid); %%%%Generate preprocessing matrix A'*A+B'*B+lambda^2*G'G
        
%         [l,u]=shangxiajie3(caiyang,grid,mx,my,h);
        if luid==1
            [l,u]=shangxiajie(caiyang,mx,my,h);
        end
      
         clear caiyang grid;
        for t=1:hasmtime
            b=quyumoni(chazhi,h)+GH;   %%%%Right end term A'*c + B'*d + lambda^2*G'*H
            xx=reshape((chazhi(2:(mx-1),2:(my-1)))',(mx-2)*(my-2),1);
            xx=tayloryuchuliCG(mx,my,b,xx,yuxishu,GG,M);  %%%%Use PCG iterations to perform calculations
            clear b;
            chazhi(2:(mx-1),2:(my-1))=(reshape(xx,(my-2),(mx-2)))';
            if luid==1
                for i=1:(mx-2)
                    for j=1:(my-2)
                        if chazhi(i+1,j+1)<l(i,j)
                            chazhi(i+1,j+1)=l(i,j);
                        elseif chazhi(i+1,j+1)>u(i,j)
                            chazhi(i+1,j+1)=u(i,j);
                        end
                    end
                end
            end
            clear xx;
        end
    end 
    jieguo=chazhi;
    clear chazhi;
end

function jieguo=taylorxishu(mx,my,GG); %%%%Generate a full coefficient matrix that can be inverted directly, this matrix is A'A + B'B +lambda^2*G'G ******
n=(mx-2)*(my-2);
e=ones(n,1);
pxia=e;
pshang=e;
mx=mx-2;
my=my-2;
pxia(my:my:(n-my))=0;
pshang((my+1):my:(n-my+1))=0;
xishu1=spdiags([pxia -2*e pshang],-1:1,n,n);%%% The matrix of coefficients of the x-direction equation A
xishu2=spdiags([e -2*e e],[-my 0 my],n,n);%%% The coefficient matrix of the equation in the y-direction B
clear pxia pshang e;
jieguo=xishu1'*xishu1+xishu2'*xishu2+GG;
clear xishu1 xishu2;

function M=duijiaoxian(mx,my); 
m=(mx-2)*(my-2);
M=12*ones(m,1);
M(1:(my-2):(m+3-my))=11;
M((my-2):(my-2):m)=11;
M(1:(my-2))=M(1:(my-2))-1;
M((m-my+3):m)=M((m-my+3):m)-1;

function x=tayloryuchuliCG(mx,my,b,x0,yuxishu,GG,M);   %%%%%Diagonal preprocessing CG iterative method yuxishu for the tridiagonal preprocessing array instead of the original equation coefficient matrices
global CGerror CGmaxtimes;
m=(mx-2)*(my-2);
mob=sqrt(sum(b.^2)/m);%%% Module of b
r=b-AP(GG,x0,mx,my,M);
clear b;
%newr=r;
ruo=sum(r.^2)/m;
x=x0;
clear x0;
k=0;
if  min(size(yuxishu))==1
    z=r./yuxishu;
else
    z=yuxishu\r;
end

while ((sqrt(ruo)>(CGerror*mob))&(k<CGmaxtimes))
    %newz=newr./M; 
    k=k+1;
    if k==1
        p=z;
    else
        beita=newrz/rz;
        p=z+beita*p;
    end
    rz=r'*z; 
    w=AP(GG,p,mx,my,M);
    alfa=rz/(p'*w);  %    w=AP(M,p,mx,my,lamda);
    x=x+alfa*p;
    r=r-alfa*w;
    if  min(size(yuxishu))==1
        z=r./yuxishu;
    else
        z=yuxishu\r;
    end
    newrz=r'*z;
    ruo=sum(r.^2)/m;
end %%%%%The iterative process requires storing r z x p w 5 matrices doing 6*m multiplications
clear r p z w yuxishu;

function pp=AP(GG,p,mx,my,M);
pp=zeros((mx-2)*(my-2),1);
pp=xiangliangweiyi2(p,mx,my,2)+xiangliangweiyi2(p,mx,my,-2)-4*xiangliangweiyi2(p,mx,my,1)-4*xiangliangweiyi2(p,mx,my,-1);
pp=pp-4*xiangliangweiyi1(p,my-2)+xiangliangweiyi1(p,2*(my-2))-4*xiangliangweiyi1(p,-(my-2))+xiangliangweiyi1(p,-2*(my-2));
pp=pp+GG*p+M.*p;

function temp=xiangliangweiyi1(p,a);
m=max(size(p));
temp=zeros(m,1);
if a>=0
    temp(1:(m-a))=p((a+1):m); %%%%top of the diagonal
elseif a<0
    temp((-a+1):m)=p(1:(m+a)); %%%%Diagonally below
end

function temp=xiangliangweiyi2(p,mx,my,a);
m=(mx-2)*(my-2);
temp=zeros(m,1);
if a==1
    p(1:(my-2):(m-my+3))=0;
elseif a==2
    p((my-1):(my-2):(m-my+3))=0;
    p(my:(my-2):(m-my+4))=0;
elseif a==(-1)
    p((my-2):(my-2):(m-my+2))=0;
elseif a==(-2)
    p((my-2):(my-2):m)=0;
    p((my-3):(my-2):(m-1))=0;
end
temp=xiangliangweiyi1(p,a);

function yuxishu=yuchulixishu(mx,my,GG,M,yuchuliid); %%%%Generating a tridiagonal array as a preprocessor, this matrix is the middle tridiagonal of A'A + B'B + lambda^2*G'G ******

if yuchuliid==1
    n=(mx-2)*(my-2);
    yuxishu=zeros(n,1);
    for i=1:n
        yuxishu(i)=M(i)+GG(i,i);
    end
else
    n=(mx-2)*(my-2);
    e=ones(n,1);
    pxia=e;
    pshang=e;
    clear e;
    mx=mx-2;
    my=my-2;
    pxia(my:my:(n-my))=0;
    pshang((my+1):my:(n-my+1))=0;
    yuxishu=spdiags([-4*pxia M -4*pshang],-1:1,n,n);
    GG=GG+yuxishu;
    clear yuxishu;
    ss=[0 (diag(GG,1))']';
    sx=[(diag(GG,-1))' 0]';
    yuxishu=spdiags([sx diag(GG) ss],-1:1,n,n); %%%% The tridiagonal of A'A+B'B+GG as a preprocessor
    clear ss sx;
end

