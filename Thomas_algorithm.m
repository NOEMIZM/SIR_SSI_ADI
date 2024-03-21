%%%% Thomas algorithm for ADI x-y, with I-2A principal diagonal and A
%%%% secundary, with Neumann
%%%% D is the solution vector

function X=Thomas_algorithm(A, D)

    n=length(D(1,:));
    I=eye(2);
    Cprime=zeros(2,2,n);
    Dprime=zeros(2,n);
    X=zeros(2,n);
    
        for i=1:n
            if i==1
                Cprime(:,:,i)=2*inv(I-2*A)*A;
                Dprime(:,i)=inv(I-2*A)*D(:,i);
            elseif i==n
                Cprime(:,:,i)=inv(I-2*A-A*Cprime(:,:,i-1))*A;
                Dprime(:,i)=inv(I-2*A-2*A*Cprime(:,:,i-1))*(D(:,i)-2*A*Dprime(:,i-1));
            else 
                Cprime(:,:,i)=inv(I-2*A-A*Cprime(:,:,i-1))*A;
                Dprime(:,i)=inv(I-2*A-A*Cprime(:,:,i-1))*(D(:,i)-A*Dprime(:,i-1));
            end
            
        end
        X(:,end)=Dprime(:,end);
        for i=length(D(1,:))-1:-1:1
            X(:,i)=Dprime(:,i)-Cprime(:,:,i)*X(:,i+1);
        end
end
