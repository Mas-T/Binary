classdef (InferiorClasses = {?sdpvar,?sym}) MDCHP < handle
% Multiple Degree Copositive Homogeneous Polynomial
%
% CHP=sum_{|alpha|=grado} z^alpha*n_alpha*coef_alpha
%      grado vector de grados del polinomio
%      n_var vector con el nº de variables combinables de forma quee
%           sum_{j=1}^{n_var(i)}z(i,j)=1
%      z monomio de n_var variables positivas
%      alpha grado de cada monomio vector dimension sum(n_var) 
%      n_alpha es la multiplicidad de z^alpha en un sumatorio completo
%      coef_alpha es el coeficiente del polinomio para el monomio z^alpha
%
%   Propiedades
%        grado  grado del polinomio
%        n_var  numero de variables
%        coef   valores de los coeficientes del polinomio  


    properties
        grado % grado del polinomio
        n_var %numero de variables
        coef %valores de los coeficientes del polinomio        
    end
    methods
        function coef_i=readcoef(thisCHP)            
            coef_i=thisCHP.coef;
        end
        function grado=readgrado(thisCHP)
            grado=thisCHP.grado;
        end
        function n_var=readn_var(thisCHP)
            n_var=thisCHP.n_var;
        end
        function newCHP=MDCHP(grado,n_var,coef)
            n_coef=1;
            for i=1:length(grado)
                n_coef=n_coef*nchoosek(grado(i)+n_var(i)-1,grado(i));
            end
            if iscell(coef)                
                if (n_coef==length(coef))
                    newCHP.coef=coef;
                    newCHP.n_var=n_var;
                    newCHP.grado=grado;
                else
                    notify(obj,'El numero de coeficientes es incorrecto');
                end
            else
                newCHP.coef=cell(1,n_coef);
                for i=1:n_coef
                    newCHP.coef{i}=coef;
                end
                 newCHP.n_var=n_var;
                 newCHP.grado=grado;
            end                
        end
        function n_coef=num_coef(thisCHP)
            % devuelve el numero de elementos del sumatorio
            %n_coef=nchoosek(thisCHP.grado+thisCHP.n_var-1,thisCHP.grado);
            %DS=CHP.multiset(thisCHP.n_var,thisCHP.grado);
            %[n_coef,~]=size(DS);
            n_coef=1;
            for i=1:length(thisCHP.grado)
                n_coef=n_coef*nchoosek(thisCHP.grado(i)+thisCHP.n_var(i)-1,thisCHP.grado(i));
            end
        end
        function DS=ConjuntoGrados(thisCHP)
            % todos los grados del polinomio
            DS=MDCHP.multiset(thisCHP.n_var,thisCHP.grado);
        end
        function coef_i=leerCoef(thisCHP,grado_monomio)
            % dado un grado devuelve su coeficiente coef_alpha
            DC=thisCHP.ConjuntoGrados();
            i=CHP.index(DC,grado_monomio);
            coef_i=thisCHP.coef{i};
        end
        function add(thisCHP,grado_monomio,coef)
            % suma el valor de coef al coeficente grado_monomio
            DC=thisCHP.ConjuntoGrados();
            i=CHP.index(DC,grado_monomio);
            thisCHP.coef{i}=thisCHP.coef{i}+coef;
        end
        function setcoef(thisCHP,grado_monomio,coef)
            % coef_grado_monomio = coef
            DC=thisCHP.ConjuntoGrados();
            i=CHP.index(DC,grado_monomio);
            thisCHP.coef{i}=coef;
        end
        function divmult(thisCHP)
            % divide cada uno de los coeficientes por su multiplicidad
            n_coef=thisCHP.num_coef();
            DS=thisCHP.ConjuntoGrados();
            for i=1:n_coef                
                thisCHP.coef{i}=thisCHP.coef{i}/CHP.permmultiset_vec(DS(i,:),thisCHP.n_var);
            end
        end
         function newCHP=minus(A,B)            
             if isa(A,'MDCHP')&& isa(B,'MDCHP') 
                 [A,B]=MDCHP.igualargrado(A,B);
                 n_coef=cellfun(@(x,y) x-y,A.coef,B.coef,'uniformOutput',false); 
                 newCHP=MDCHP(A.grado,A.n_var,n_coef);
             elseif isa(A,'MDCHP')
                 n_coef=cellfun(@(x) x-B,A.coef,'uniformOutput',false);
                 newCHP=MDCHP(A.grado,A.n_var,n_coef);
             else
                 n_coef=cellfun(@(x) A-x,B.coef,'uniformOutput',false);
                 newCHP=MDCHP(B.grado,B.n_var,n_coef);
             end           
        end
        function newCHP=uminus(A)            
             n_coef=cellfun(@(x) -x,A.coef,'uniformOutput',false);
             newCHP=MDCHP(A.grado,A.n_var,n_coef);             
        end
        function newCHP=uplus(A)            
             n_coef=cellfun(@(x) +x,A.coef,'uniformOutput',false);
             newCHP=MDCHP(A.grado,A.n_var,n_coef);             
        end
        function newCHP=mrdivide(A,B)
             if isa(A,'MDCHP')&& ~isa(B,'MDCHP')
                n_coef=cellfun(@(x) x/B,A.coef,'uniformOutput',false);
                newCHP=MDCHP(A.grado,A.n_var,n_coef);             
             else
                error("Undefined operator '/' for input arguments of type 'MDCHP'");
             end
        end
        function Polya(thisCHP,nuevo_grado)
            % realiza la expansion de polya hasta el grado nuevo_grado
            q=nuevo_grado;
            r=thisCHP.n_var;
            d=thisCHP.grado;
            Q=thisCHP.coef;
            MS_Q=thisCHP.ConjuntoGrados();
            MS=MDCHP.multiset(r,q);
            [n,~]=size(MS);
            tilde_Q=cell(1,n);
            for i=1:n
                tilde_Q{i}=0;
                gamma=MS(i,:);
                SMS=MDCHP.submultiset_vec(gamma,d,r);
                [sn,~]=size(SMS);
                for j=1:sn
                    alpha=SMS(j,:);
                    tilde_Q{i}=tilde_Q{i}+MDCHP.permmultiset_vec(alpha,r)*MDCHP.permmultiset_vec((gamma-alpha),r)*Q{MDCHP.index(MS_Q,alpha)};
                end
                tilde_Q{i}=tilde_Q{i}/MDCHP.permmultiset_vec(gamma,r);
            end
            thisCHP.coef=tilde_Q;
            thisCHP.grado=q;
        end
        function newCHP=mtimes(A,B) 
            if isa(A,'MDCHP')&& isa(B,'MDCHP')               
                r=A.n_var;
                d=A.grado;
                q=d+B.grado;
                Q=A.coef;
                L=B.coef;
                MS_Q=A.ConjuntoGrados();
                MS_L=B.ConjuntoGrados();
                MS=MDCHP.multiset(r,q);
                [n,~]=size(MS);
                tilde_Q=cell(1,n);
                for i=1:n
                    tilde_Q{i}=0;
                    gamma=MS(i,:);
                    SMS=MDCHP.submultiset_vec(gamma,d,r);
                    [sn,~]=size(SMS);
                    for j=1:sn
                        alpha=SMS(j,:);
                        tilde_Q{i}=tilde_Q{i}+MDCHP.permmultiset_vec(alpha,r)*MDCHP.permmultiset_vec(gamma-alpha,r)*Q{MDCHP.index(MS_Q,alpha)}*L{MDCHP.index(MS_L,gamma-alpha)};
                    end
                    tilde_Q{i}=tilde_Q{i}/MDCHP.permmultiset_vec(gamma,r);
                end
                newCHP=MDCHP(q,r,tilde_Q);
            elseif isa(A,'MDCHP')
                n_coef=cellfun(@(x) x*B,A.coef,'uniformOutput',false);
                newCHP=MDCHP(A.grado,A.n_var,n_coef);                
            else
                n_coef=cellfun(@(x) A*x,B.coef,'uniformOutput',false);
                newCHP=MDCHP(B.grado,B.n_var,n_coef); 
            end
                    
        end
        function newCHP=plus(A,B)            
             if isa(A,'MDCHP')&& isa(B,'MDCHP') 
                 [A,B]=MDCHP.igualargrado(A,B);
                 n_coef=cellfun(@(x,y) x+y,A.coef,B.coef,'uniformOutput',false); 
                 newCHP=MDCHP(A.grado,A.n_var,n_coef);
             elseif isa(A,'MDCHP')
                 n_coef=cellfun(@(x) x+B,A.coef,'uniformOutput',false);
                 newCHP=MDCHP(A.grado,A.n_var,n_coef);
             else
                 n_coef=cellfun(@(x) A+x,B.coef,'uniformOutput',false);
                 newCHP=MDCHP(B.grado,B.n_var,n_coef);
             end           
        end
        function newCHP=double(A)
             n_coef=cellfun(@(x) double(x),A.coef,'uniformOutput',false);
             newCHP=MDCHP(A.grado,A.n_var,n_coef);
        end
        function newCHP=vertcat(varargin)
            Aux=varargin{1};
            i=1;
            while  ~isa(Aux,'MDCHP')
                i=i+1;
                Aux=varargin{i};
            end
            nvar=Aux.n_var;
            grad=Aux.grado;
            Aux=varargin{1};
            if isa(Aux,'MDCHP')
                Prev=MDCHP(Aux.grado,Aux.n_var,Aux.coef);
            else
                Prev=MDCHP(grad,nvar,Aux);
            end
            for i=2:nargin
                Aux=varargin{i};
                if isa(Aux,'MDCHP')
                    Act=MDCHP(Aux.grado,Aux.n_var,Aux.coef);
                else
                    Act=MDCHP(grad,nvar,Aux);
                end
                [Prev,Act]=MDCHP.igualargrado(Prev,Act);
                n_coef=cellfun(@(x,y) [x; y],Prev.coef,Act.coef,'uniformOutput',false); 
                Prev=MDCHP(Act.grado,nvar,n_coef);
            end
            newCHP=Prev;
        end
        function newCHP=horzcat(varargin)
            Aux=varargin{1};
            i=1;
            while  ~isa(Aux,'MDCHP')
                i=i+1;
                Aux=varargin{i};
            end
            nvar=Aux.n_var;
            grad=Aux.grado;
            Aux=varargin{1};
            if isa(Aux,'MDCHP')
                Prev=MDCHP(Aux.grado,Aux.n_var,Aux.coef);
            else
                Prev=MDCHP(grad,nvar,Aux);
            end
            for i=2:nargin
                Aux=varargin{i};
                if isa(Aux,'MDCHP')
                    Act=MDCHP(Aux.grado,Aux.n_var,Aux.coef);
                else
                    Act=MDCHP(grad,nvar,Aux);
                end
                [Prev,Act]=MDCHP.igualargrado(Prev,Act);
                n_coef=cellfun(@(x,y) [x y],Prev.coef,Act.coef,'uniformOutput',false); 
                Prev=MDCHP(Act.grado,nvar,n_coef);
            end
            newCHP=Prev;
        end
        function newCHP=ctranspose(A)
            n_coef=cellfun(@(x) x',A.coef,'uniformOutput',false);
            newCHP=MDCHP(A.grado,A.n_var,n_coef);
        end
        function newCHP=le(A,B)
            if isa(A,'MDCHP')&& isa(B,'MDCHP')
                [A,B]=MDCHP.igualargrado(A,B);
                n_coef=cellfun(@(x,y) x<=y,A.coef,B.coef,'uniformOutput',false);
                newCHP=MDCHP(A.grado,A.n_var,n_coef);
            elseif isa(A,'MDCHP')
                n_coef=cellfun(@(x) x<=B,A.coef,'uniformOutput',false);
                newCHP=MDCHP(A.grado,A.n_var,n_coef);
            else
                n_coef=cellfun(@(x) A<=x,B.coef,'uniformOutput',false);
                newCHP=MDCHP(B.grado,B.n_var,n_coef);
            end
        end
        function newCHP=ge(A,B)
            if isa(A,'MDCHP')&& isa(B,'MDCHP')
                [A,B]=MDCHP.igualargrado(A,B);
                n_coef=cellfun(@(x,y) x>=y,A.coef,B.coef,'uniformOutput',false);
                newCHP=MDCHP(A.grado,A.n_var,n_coef);
            elseif isa(A,'MDCHP')
                n_coef=cellfun(@(x) x>=B,A.coef,'uniformOutput',false);
                newCHP=MDCHP(A.grado,A.n_var,n_coef);
            else
                n_coef=cellfun(@(x) A>=x,B.coef,'uniformOutput',false);
                newCHP=MDCHP(B.grado,B.n_var,n_coef);
            end
        end
        function coef=array(thisCHP)
            %devuelve los coeficientes en un array
            coef=[];
            for i=1:thisCHP.num_coef()
                coef=[coef; thisCHP.coef{i}];
            end
        end
        function b=isnan(thisCHP)
            b=isnan(thisCHP.coef{1});
        end
        function [n,m]=size(thisCHP)
            [n,m]=size(thisCHP.coef{1});
        end
        function y=eval(thisCHP,x)
            y=0;
            r=thisCHP.n_var;
            DS=thisCHP.ConjuntoGrados();
            for i=1:size(DS,1)
                y=y+prod(x.^DS(i,:))*thisCHP.coef{i}*MDCHP.permmultiset_vec(DS(i,:),r);
            end
        end
    end
    methods (Static)  
        function [An,Bn]=igualargrado(A,B)
            if A.grado==B.grado
                An=A;
                Bn=B;
            else
                n_grado=max(A.grado,B.grado);
                Bn=MDCHP(B.grado,B.n_var,B.coef);
                Bn.Polya(n_grado);
                An=MDCHP(A.grado,A.n_var,A.coef);
                An.Polya(n_grado);
            end
        end
        function n=permmultiset(alpha)
            n=1;
            for i=1:length(alpha)
                n=n*factorial(sum(alpha{i}));
                for j=1:length(alpha{i})
                    n=n/factorial(alpha{i}(j));
                end
            end
        end
        function n=permmultiset_vec(alpha,r)
            n=1;
            ind1=1;            
            for i=1:length(r)
                ind2=ind1+r(i)-1;
                aux=alpha(ind1:ind2);
                ind1=ind2+1;
                n=n*factorial(sum(aux));
                for j=1:length(aux)
                    n=n/factorial(aux(j));
                end
            end
        end
        function MS=multiset(r,q)
            MS=MDCHP.smultiset(r(1),q(1));
            for i=2:length(r)
                SMS=MDCHP.smultiset(r(i),q(i));
                [n,~]=size(SMS);
                [m,~]=size(MS);
                MS=kron(MS,ones(n,1));
                SMS=kron(ones(m,1),SMS);
                MS=[MS SMS];
            end
        end
        function MS=smultiset(r,q)
            % r es el numero de reglas o variables
            % q es la multiplicidad o grado de los monomios del polinomio      
            if (q==0)
                MS=zeros(1,r);
            else
                MS=[];
                if (r>1)
                    for i=q:-1:0
                        MSaux=MDCHP.smultiset(r-1,q-i);
                        [n,~]=size(MSaux);
                        MSfila=[i*ones(n,1) MSaux];
                        MS=[MS; MSfila];
                    end
                else
                    MS=q;
                end
                
            end
        end
        function i=index(Multiset,alpha)
            % indice de alpha en el multiset
            i=1;
            if iscell(alpha)
                alpha=cell2mat(alpha);
            end
            [n,~]=size(Multiset);
            while sum(Multiset(i,:)~=alpha)
                if (i>n)
                    break
                end
                i=i+1;
            end        
        end
        function SMS=submultiset(gamma,d)
            % gamma es el indice del que extraemos el subconjunto
            % d es el orden del monomio extraido
            r=zeros(size(gamma));
            for i=1:length(gamma)
                r(i)=length(gamma{i});
            end
            MS=MDCHP.multiset(r,d);
            [n,~]=size(MS);
            SMS=[];
            for i=1:n
                if (MS(i,:)<=cell2mat(gamma))
                    SMS=[SMS;MS(i,:)];
                end
            end
        end
        function SMS=submultiset_vec(gamma,d,r)
            % gamma es el indice del que extraemos el subconjunto
            % d es el orden del monomio extraido
            
            MS=MDCHP.multiset(r,d);
            [n,~]=size(MS);
            SMS=[];
            for i=1:n
                if (MS(i,:)<=gamma)
                    SMS=[SMS;MS(i,:)];
                end
            end
        end
        function n_coef=number_coef(grado,n_var)
            % devuelve el numero de elementos del sumatorio
            %n_coef=nchoosek(thisCHP.grado+thisCHP.n_var-1,thisCHP.grado);
            %DS=CHP.multiset(thisCHP.n_var,thisCHP.grado);
            %[n_coef,~]=size(DS);
            n_coef=1;
            for i=1:length(grado)
                n_coef=n_coef*nchoosek(grado(i)+n_var(i)-1,grado(i));
            end
        end
    end
end