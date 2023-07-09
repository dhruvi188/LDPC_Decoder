clear all;
m = 3792;
n = 5056;

H1=load("Hmatrix.mat");
H=H1.H;

dc=0;
    for i=1:n
        if(H(1,i)==1)
            dc=dc+1;
        end
    end %found dc
    
dv=0;
    for i=1:m
        if(H(i,1)==1)
            dv=dv+1;
        end
    end %found dv
    
    V_dv=zeros(n, dv); % ith index stores CNs connected to VN i
    C_dc=zeros(m,dc); %ith index tores VNs connected to CN i
    
    for i=1:m %preparing tanner graph(dc)
        k=1;
        for j=1:n
            if(H(i,j)==1)
                C_dc(i,k)=j;
                k=k+1;
            end
        end
    end
    
    for i=1:n %preparing tanner graph(dv)
        k=1;
        for j=1:m
            if(H(j,i)==1)
                V_dv(i,k)=j;
                k=k+1;
            end
        end
    end

prompt="input probability of error";
po = input(prompt);
success=0;
errorProbab = zeros(51,1);
errorProbab(1) = po;
for simIdx = 1:100

vn = -1*(rand(n,1)>(1-po));


vnprob = ones(n,1);
for i=1:n
    if(vn(i)==0)
        vnprob(i)=0;
    elseif(vn(i)==1)
        vnprob(i)=1;
    elseif(vn(i)==-1)
        vnprob(i)=0.5;
    end
end

H_prob = H;

for i=1:n
    for j=1:dv
            if(V_dv(i,j)>0)
                H_prob(V_dv(i,j),i) = vnprob(i);
            end
    end
end

eraser = ones(50,1);




for t=1:50
    
    num=0;
    for i=1:n
        if(vn(i)==-1)
               num=num+1;
        end
    end
    errorProbab(t+1) = errorProbab(t+1) + num/(n*100);
    vncopy=vn;
    for i=1:m
        vntemp = zeros(dc,1);
        for j=1:dc
            temp=1;
            for k=1:dc
                if(k~=j)
                    temp=temp*(1-2*(H_prob(i,C_dc(i,k))));
                end
            end
            vntemp(j)=0.5 +0.5*temp;
        end
        for it = 1:dc
            H_prob(i,C_dc(i,it))=vntemp(it);
        end
    end
    
    for i=1:n
        vntemp=zeros(dv,1);
        for j=1:dv
            p1=1;
            p0=1;
            for k=1:dv
                if(k~=j)
                    p1=p1*(1-H_prob(V_dv(i,k),i));
                    p0=p0*(H_prob(V_dv(i,k),i));
                end
            end
            p1 = p1*vnprob(i);
            p0 = p0*(1-vnprob(i));
            vntemp(j) = p1/(p1+p0);
        end
        for k=1:dv
            H_prob(V_dv(i,k), i) = vntemp(k);
        end
    end

    for i=1:n
        transmitted1=1;
        transmitted0=1;
        for j=1:dv
            transmitted1 = transmitted1*(H_prob(V_dv(i,j),i));
            transmitted0 = transmitted0*(1-H_prob(V_dv(i,j),i));
        end
        transmitted1 = transmitted1*vnprob(i);
        transmitted0 = transmitted0*(1-vnprob(i));
        if(transmitted0>transmitted1)
            vnprob(i)=0;
            vn(i)=0;
        elseif(transmitted1>transmitted0)
            vnprob(i)=1;
            vn(i) =1;
        end
    end
    %if(vncopy==vn)
     %   break;
    %end
    %if(~(any(vn)))
     %   break;
    %end
    
end
if(~any(vn))
    success=success+1;
end

end

idx=zeros(51,1);
for i=1:51
    idx(i)=i;
end
hold on;
plot(idx, errorProbab);

analytical = zeros(51,1);
analytical(1) = po;
for i=1:50
    analytical(i+1) = po*(1-((1-analytical(i))^(dc-1)))^(dv-1);
end

sucessRatio = success/100;

plot(idx,analytical)
grid on;
legend('decoding algorithm convergence','analytical convergence');
hold off;
