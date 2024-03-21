%%% X is a 2 x n x n matrix
function R=SIR_Reaction(X, beta, gamma)
    n=length(X(1,:,1));
    U=X(1,:,:);
    V=X(2,:,:);
    R1=-beta*U.*V;
    R2=beta*U.*V-gamma*V;
    R=zeros(2,n,n);
    R(1,:,:)=R1;
    R(2,:,:)=R2;
end