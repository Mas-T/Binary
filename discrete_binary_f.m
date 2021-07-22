
function [PA,PB,alpha2]=discrete_binary(A,B,T)
nu=size(B,2);
r=rank(A);
[V,D]=eig(A); %diagonalize A
if ~isreal(D)
    [V,D] = cdf2rdf(V,D); %if D is not real, transforms the current diagonal 
%matrix into another one with the complex-conjugate values turned into real 2x2 blocks
end
n=length(D);
d=zeros(1,n);
alpha2=cell(r,1);
vars=2*ones(1,r);
k=1; %counter for each fuzzy model
PB=zeros(n,nu);
j=1; %counter for each diagonal value
PD=zeros(n,n);
%Now, it iterates trhu all D values:
while j<n
        if (D(j,j+1)==0)&&(D(j,j)~=0) %if that value is real:
           d(j)=(1-exp(-D(j,j)*T))/D(j,j)/T; 
           dbk1=blkdiag(zeros(j-1,j-1),1,zeros(n-j,n-j)); %lower limit 
           dbk2=blkdiag(zeros(j-1,j-1),d(j),zeros(n-j,n-j)); %upper limit
           alpha2{k}=@(delta) ((1-exp(-D(j,j)*T*delta))/D(j,j)./T/delta-1)/(d(j)-1); %calculates the membership functions that are delta dependant.
           PD=PD+MDCHP([zeros(1,k-1) 1 zeros(1,r-k)],vars, {dbk1 dbk2}); %adds the previous info into a block to the corresponding position in PD
           k=k+1;
           j=j+1;
        elseif (D(j,j+1)==0)&&(D(j,j)==0)%if that value is 0, a pole on the origin.
            PD=PD+blkdiag(zeros(j-1,j-1),1,zeros(n-j,n-j)); %creates a new position in PD but adds no info to it, since there's not any fuzzy system.
            j=j+1;
        else %if that value is complez:
            a11=D(j,j);
            a12=D(j,j+1);
            aQ=a11^2+a12^2;
            [f_min,f_max,alpha2{k}]=fuzzy_complex_poles_f_1(a11,a12,T); %extracts the functions f(x) & f(x+1) to operate the limits. 
            %Also extracts the membership functions for that value.
            d1=[-a11/aQ a12/aQ;-a12/aQ -a11/aQ]*f_min; %lower limit calculation
            d2=[-a11/aQ a12/aQ;-a12/aQ -a11/aQ]*f_max;%upper limit calculation
            dbk1=blkdiag(zeros(j-1,j-1),d1,zeros(n-j-1,n-j-1)); %lower limit on a matrix form
            dbk2=blkdiag(zeros(j-1,j-1),d2,zeros(n-j-1,n-j-1)); %upper limit on a matrix form
            PD=PD+MDCHP([zeros(1,k-1) 1 zeros(1,r-k)],vars, {dbk1 dbk2}); %adds the previous info into a block to the corresponding position in PD
            k=k+1;
            [f_min,f_max,alpha1]=fuzzy_complex_poles_f_2(a11,a12,T); %now calculates the same for the complex-conjugate
            d1=[a12/aQ a11/aQ;-a11/aQ a12/aQ]*f_min;
            d2=[a12/aQ a11/aQ;-a11/aQ a12/aQ]*f_max;
            alpha2{k}=@(x) alpha1(x);            
            dbk1=blkdiag(zeros(j-1,j-1),d1,zeros(n-j-1,n-j-1));
            dbk2=blkdiag(zeros(j-1,j-1),d2,zeros(n-j-1,n-j-1));
            j=j+1;
            PD=PD+MDCHP([zeros(1,k-1) 1 zeros(1,r-k)],vars, {dbk1 dbk2});%and adds it to PD
            k=k+1;
            j=j+1;
        end
end
if j==n; %This section only applyes for the last item of the diagonal since the counter may left it unsolved under certain conditions 
    %and does exactly the same as shown above.
    if (D(n,n)~=0)%If the value is real
           d(j)=(1-exp(-D(j,j)*T))/D(j,j)/T; 
           dbk1=blkdiag(zeros(j-1,j-1),1,zeros(n-j,n-j));
           dbk2=blkdiag(zeros(j-1,j-1),d(j),zeros(n-j,n-j));
           alpha2{k}=@(delta) ((1-exp(-D(j,j)*T*delta))/D(j,j)./T/delta-1)/(d(j)-1);
           PD=PD+MDCHP([zeros(1,k-1) 1 zeros(1,r-k)],vars, {dbk1 dbk2});
           k=k+1;
    else % or if that value is a zero.
           PD=PD+blkdiag(zeros(j-1,j-1),1,zeros(n-j,n-j));
        end        
end 

%Now, with PD defined, the discrete matrix PB can be solved.
coef=PD.coef;
unos=ones(1,r);
vars2=2*ones(1,nu*r);
PB=zeros(n,nu);
aux=alpha2;
alpha2={};
for i=1:nu %for each control action, a different PB is required
    PDi=MDCHP([zeros(1,r*(i-1)) unos zeros(1,r*(nu-i))],vars2,coef); %Creates an auxiliar array matrix where each position is 
    %an info block that will contain the PB fuzzy information for each control action.
    Bi=[zeros(n,i-1) B(:,i) zeros(n,nu-i)]; %Trims B into each control action(i) and completes with zeros the rest of the matrix
    PB=PB+T*expm(A*T)*V*PDi*inv(V)*Bi; %Calculates the corresponding PB for that specific control action and adds it to the complete PB matrix.
    alpha2=[alpha2; aux]; %Adds the alphas for that control action on a separate line, creating the alphas matrix where each line goes for each PB(i)
end
PA=expm(A*T); %PA is straigthfoward.
end