%% This file will generate the guess of the iSINDYc left hand side.
% Last Updated: 2019/04/21
% Coded By: K
function [Data,Sym_Struct]=GuessLib(X,dX,iter,u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order, Integral, dt)
%% First get the size of the X matrix, determin the data length and the number of variables we have.
[Data_Length,Variable_Number]=size(X);
[~,Variable_Number_dX]=size(dX);
[~,Variable_Number_u]=size(u);
% Also create the symbolic variable
Symbol=sym('z',[Variable_Number,1]);
Symbol_dX=sym('dz',[Variable_Number,1]);
Symbol_u=sym('u',[Variable_Number_u,1]);

%% Now according the Highest Polynomial Order entered, we will calculate the data matrix.
Data=[];
Index=1;

% Form basis vector
Basis=[X(:,1) X(:,2) X(:,1)-X(:,2) X(:,1)-2*X(:,2) 2*X(:,1)-X(:,2) 2*X(:,1)-2*X(:,2)];
Basis_Sym=[Symbol(1) Symbol(2) Symbol(1)-Symbol(2) Symbol(1)-2*Symbol(2) 2*Symbol(1)-Symbol(2) 2*Symbol(1)-2*Symbol(2)];

% Add dx term, special cases during iSINDy, index 1 and 2. Same on both RHS
% and LHS
Data(:,Index)=dX(:,1);
Sym_Struct{1,Index}=Symbol_dX(iter);
Index=Index+1;

Data(:,Index)=dX(:,1).*cos(Basis(:,2)).^2;
Sym_Struct{1,Index}=Symbol_dX(iter)*cos(Basis_Sym(1,2)).^2;
Index=Index+1;

%old
% for i=3:size(Basis,2)
%    Data(:,Index)=dX(:,1).*cos(Basis(:,i)).^2;
%    Sym_Struct{1,Index}=Symbol_dX(iter)*cos(Basis_Sym(1,i)).^2;
%    Index=Index+1;
% end

%not dx terms

Data(:,Index)=ones(Data_Length,1);
Sym_Struct{1,Index}=1;

%Could probably be improved, with integration methods more specific to the
%library terms
if Integral == 1
    %prepare compensatory term, only used for iter > 2
    comp = -sin(2*Basis(:,2)).*X(:,4).*X(:,iter);

    %integrate all first order or lower derivatives and comp, trapezoidal TODO.
    for i = 2:length(Data(:,1))
        %create temporary integal value lists
        % TID = Data;
        % TIC = comp;
        % TID(i,:) = TID(i-1,:) + (Data(i-1,:) + Data(i,:))/2.*dt;
        % TIC(i,:) = TIC(i-1,:) + (comp(i-1,:) + comp(i,:))/2.*dt;
        % 
        % 
        % Data = TID;
        % comp = TIC;

        Data(i,:) = Data(i-1,:) + Data(i,:).*dt;
        comp(i,:) = comp(i-1,:) + comp(i,:).*dt;
    end

    %integrate the second order ones
    if iter > 2
        %first inedx fix
        Data(:,1) = X(:,iter);
        %second index fix
        Data(:,2) = X(:,iter).*cos(Basis(:,2)).^2 - comp;
    end

end
end

