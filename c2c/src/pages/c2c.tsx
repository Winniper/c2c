import Card from '../components/Card'

const C2C = () => {
  const bigM = `\
format short
clear all
clc

% Define the cost vector (objective function coefficients)
% Artificial variables are penalized with -10000
Cost = [-2 -1 0 0 -10000 -10000 0]; % Extra column for RHS

% Constraint matrix A with slack and artificial variables included
A = [3  1  0  0  1  0  3;
     4  3 -1  0  0  1  6;
     1  2  0  1  0  0  3];

% Basic variable indices (artificial and slack variables in basis)
BV = [5 6 4];  % Indices correspond to A1, A2, S2

% Initial Zj - Cj row
ZjCj = Cost(BV) * A - Cost;

% Display initial tableau
zcj = [Cost; ZjCj; A];
bigmtable = array2table(zcj);
bigmtable.Properties.VariableNames(1:size(zcj,2)) = {'x_1','x_2','s_1','s_2','A_1','A_2','sol'};
disp('Initial Simplex Tableau:')
disp(bigmtable)

% Simplex iterations
RUN = true;
while RUN
    ZC = ZjCj(1:end-1); % Exclude solution column

    if any(ZC < 0)
        fprintf('\nThe current BFS is not optimal\n');

        % Entering variable (most negative Zj-Cj)
        [~, pvt_col] = min(ZC);
        fprintf('Entering Column = %d\n', pvt_col);

        sol = A(:, end);         % Solution column
        Column = A(:, pvt_col);  % Pivot column

        % Check for unboundedness
        if all(Column <= 0)
            error('LPP is unbounded');
        else
            % Compute ratios for the minimum ratio test
            ratio = zeros(size(A,1), 1);
            for i = 1:size(A,1)
                if Column(i) > 0
                    ratio(i) = sol(i) / Column(i);
                else
                    ratio(i) = inf;
                end
            end

            % Determine leaving variable
            [~, pvt_row] = min(ratio);
            fprintf('Leaving Row = %d\n', pvt_row);
        end

        % Pivot operations
        BV(pvt_row) = pvt_col;  % Update basis
        pvt_key = A(pvt_row, pvt_col);
        A(pvt_row, :) = A(pvt_row, :) / pvt_key;

        for i = 1:size(A,1)
            if i ~= pvt_row
                A(i, :) = A(i, :) - A(i, pvt_col) * A(pvt_row, :);
            end
        end

        % Update Zj - Cj row
        ZjCj = Cost(BV) * A - Cost;

        % Display updated tableau
        ZCj = [Cost; ZjCj; A];
        TABLE = array2table(ZCj);
        TABLE.Properties.VariableNames(1:size(ZCj,2)) = {'x_1','x_2','s_1','s_2','A_1','A_2','sol'};
        disp('Updated Simplex Tableau:')
        disp(TABLE)

    else
        RUN = false;
        fprintf('\nCurrent BFS is Optimal.\n');
        Z_optimal = Cost(BV) * A(:, end);
        fprintf('Optimal Value of Z = %f\n', Z_optimal);
    end
end
`;

const simplex =`\
% Max. z = x1 + 2x2, subject to − x1 + x2 ≤ 1, x1 + x2 ≤ 2, x1, x2 ≥ 0.

% Max. z = 4x1 + 6x2 + 3x3 + x4
% subject to
% x1 + 4x2 + 8x3 + 6x4 ≤ 11, 4x1 + x2 + 2x3 + x4 ≤ 7, 2x1 + 3x2 + x3 + 2x4 ≤ 2, x1, x2, x3 ≥ 0.
% (Ans: x1 = 1/3, x3 = 4/3s2 = 3; z = 16/3)

% Min. z = −3/4x4 + 20x5 − 1/2x6 + 6x7 
% subject to
% x1 + 1/4x4 − 8x5 − x6 + 9x7 = 0, x2 + 1/2x4 − 12x5 − 1/6x6 + 3x7 = 0, x3 + x6 = 1
% (Ans: x1 = 3/4, x4 = 1, x6 = 1; z = −5/4)

clc;
clear all;
format short;

% Number of decision variables
N = 4;

% Objective Function Coefficients (Negative for maximization)
C = [4 6 3 1]; 

% Constraint Coefficients
A = [1 4 8 6;4 1 2 1;2 3 1 2];

% Right-hand side of constraints
b = [11;7;2];

% Slack Variables (Identity matrix)
s = eye(3);

% Constructing the Simplex Table
A = [A s b];

% Cost row initialization
cost = zeros(1, size(A,2));
cost(1:N) = C;

% Basic Variables (Initially Slack Variables)
bv = [5 6 7]; 

% Compute Zj - Cj Row
zjcj = cost(bv) * A - cost;

% Initial Simplex Table

simptable = array2table([zjcj; A])
simptable.Properties.VariableNames(1:size(zjcj,2)) = {'x_1','x_2','x_3','x_4','s_1','s_2','s_3','sol'}
RUN= true
while RUN
if any(zjcj<0)
    zc=zjcj(1:end-1)% from x1 upto s3 leaving last col soln
 [Enter,pivot_col]=min(zc) % entering pivot col
 if all(A(:,pivot_col)<=0)
    error('LPP is Unbounded all enteries are <=0 in column %d',pivot_col);
   else
 sol=A(:,end)
 column=A(:,pivot_col)
 for i=1:size(column,1)
     if column(i)>0
         ratio(i)=sol(i)/column(i);
     else
         ratio(i)=inf;
     
     end
 end
     [leaving_val, pivot_row]=min(ratio)
 end
bv(pivot_row)=pivot_col
pivot_key=A(pivot_row, pivot_col)
A(pivot_row,:)=A(pivot_row,:)./pivot_key
for i=1:size(A,1)
    if i~=pivot_row
        A(i,:)=A(i,:)-A(i, pivot_col).*A(pivot_row,:)
    end
end
zjcj=zjcj-zjcj(pivot_col).*A(pivot_row,:)
    zcj=[zjcj;A]
    table=array2table(zcj)
    table.Properties.VariableNames(1:size(zcj,2))= {'x_1','x_2','x_3','x_4','s_1','s_2','s_3','sol'}
  
   bfs=zeros(1,size(A,2))
bfs(bv)= A(:,end)
bfs(end)=sum(bfs .* cost)
current_bfs=array2table(bfs)
current_bfs.Properties.VariableNames(1:size(current_bfs,2))={'x_1','x_2','x_3','x_4','s_1','s_2','s_3','sol'}
    else
    RUN=false;
    fprintf('The current BFS is optimal \n')
end
end
`;

  const dualSimplex = `\
clc;
clear all;
format short;

% Number of decision variables
N = 2;

% Objective Function Coefficients (Negative for maximization)
C = [-5 -6]; 

% Constraint Coefficients
A = [-1 -1 ;-4 -1];

% Right-hand side of constraints
b = [-2;-4];

% Slack Variables (Identity matrix)
s = eye(2);

% Constructing the Simplex Table
A = [A s b];

% Cost row initialization
cost = zeros(1, size(A,2));
cost(1:N) = C;

% Basic Variables (Initially Slack Variables)
bv = [3 4]; 

% Compute Zj - Cj Row
zjcj = cost(bv) * A - cost;

% Initial Simplex Table

simptable = array2table([zjcj; A])
simptable.Properties.VariableNames(1:size(zjcj,2)) = {'x_1','x_2','s_1','s_2','sol'}
RUN= true
while RUN
    if any(A(:,end) < 0) % Check infeasibility
        fprintf('Current BFS is not feasible\n');
        [minVal, pivot_row] = min(A(:,end)); % Most negative solution
        fprintf('Leaving row: %d\n', pivot_row);
        
        row = A(pivot_row, 1:end-1);
        z_row = zjcj(1:end-1);
        
        ratio = inf(1, length(row));
        for j = 1:length(row)
            if row(j) < 0
                ratio(j) = z_row(j) / row(j);
            end
        end
        
        [minRatio, pivot_col] = min(abs(ratio)); % min ratio
        fprintf('Entering column: %d\n', pivot_col);
        
        bv(pivot_row) = pivot_col;
        pivot_key = A(pivot_row, pivot_col);
        A(pivot_row,:) = A(pivot_row,:) ./ pivot_key;
        
        for i = 1:size(A,1)
            if i ~= pivot_row
                A(i,:) = A(i,:) - A(i,pivot_col) .* A(pivot_row,:);
            end
        end
        
        zjcj = zjcj - zjcj(pivot_col) .* A(pivot_row,:);
        zcj = [zjcj; A];
        table = array2table(zcj);
        table.Properties.VariableNames(1:size(zcj,2))={'x_1','x_2','s_1','s_2','sol'};
        
        bfs = zeros(1, size(A,2));
        bfs(bv) = A(:,end);
        bfs(end) = sum(bfs .* cost);
        current_bfs = array2table(bfs);
        current_bfs.Properties.VariableNames(1:size(current_bfs,2)) = {'x_1','x_2','s_1','s_2','sol'}
        
    else
        RUN = false;
        fprintf('The current BFS is Optimal & Feasible \n');
    end
end
`;

const LCM = `\
clc
clear all
format short
% Matlab Code of Least Cost Method (LCM)
% Input Information
%% Input Phase
Cost=[2 10 4 5 ;6 12 8 11 ;3 9 5 7 ]
A=[12 25 20]
B=[25 10 15 5]

%% To check unbalanced/balanced Problem
if sum(A)==sum(B)
    fprintf('Given Transportation Problem is Balanced \n')
else
   fprintf('Given Transportation Problem is Unbalanced \n') 
   if sum(A)<sum(B)
       Cost(end+1,:)=zeros(1,size(B,2))
       A(end+1)=sum(B)-sum(A)
   elseif sum(B)<sum(A)
   Cost(:,end+1)=zeros(1,size(A,2))
       B(end+1)=sum(A)-sum(B)  
   end
end

ICost=Cost
X=zeros(size(Cost))   % Initialize allocation
[m,n]=size(Cost)      % Finding No. of rows and columns
BFS=m+n-1             % Total No. of BFS

%% Finding the cell(with minimum cost) for the allocations
for i=1:size(Cost,1)  %First loop goes row by row
    for j=1:size(Cost,2)  %Second loop explores every column for each row
hh=min(Cost(:))   %Find the smallest cost in the whole matrix (flattened into a single column)
[Row, Col]=find(hh==Cost) %Find where that cheapest cell is in the matrix."
x11=min(A(Row),B(Col))  % What’s the maximum we can ship here without exceeding supply or demand
[Value,index]=max(x11)        %Among all those minimums, find the biggest shipment you can make.
ii=Row(index)       % Store selected row
jj=Col(index)        % Store selected column 
y11=min(A(ii),B(jj))        % Allocate as much as possible to this cell
X(ii,jj)=y11           % Put this value in the allocation matrix.
A(ii)=A(ii)-y11         % Reduce supply after allocation
B(jj)=B(jj)-y11         % Reduce demand after allocation
Cost(ii,jj)=inf        %Set used cell cost to infinity so it's skipped later.
    end
end

%% Print the initial BFS
fprintf('Initial BFS =\n')
IBFS=array2table(X)
disp(IBFS)

%% Check for Degenerate and Non Degenerate
TotalBFS=length(nonzeros(X))
if TotalBFS==BFS
    fprintf('Initial BFS is Non-Degenerate \n')
else
    fprintf('Initial BFS is Degenerate \n')
end


%% Compute the Initial Transportation cost
InitialCost=sum(sum(ICost.*X))
fprintf('Initial BFS Cost is = %d \n',InitialCost)
`;

const SD = `\
syms x y
f1= x^2+ y^2
fx= inline(f1)
fobj= @(x) fx(x(:,1), x(:,2));

g1=gradient(f1)
gx=inline(g1)
gobj= @(x) gx(x(:,1), x(:,2));

h1=hessian(f1)
hx= inline(h1)

x0=[1 2]
tol=0.05
maxiter=2
iter=0

x=[]

while norm(gobj(x0))> tol && iter< maxiter

    x=[x; x0]
    s=-gobj(x0);
    h=hx(x0);
    lambda= s' * s/ (s'* h*s)
    xnew= x0+ lambda* s'
    x0=xnew
    iter=iter+1
end
fprintf('Optimal Solution x = [%f, %f]\n', x0(1), x0(2));
fprintf('Optimal value f(x) = %f \n', fobj(x0));
`;

  const copyToClipboard = (text: string) => {
    navigator.clipboard.writeText(text.replace(/\n/g, "\r\n"))
  };

  return (
    <div className='bg-[#0D1117] h-dvh w-auto flex justify-center items-center'>
      <div className='flex gap-5 pt-24 pb-10'>
        <Card classname="bg-[#584F44] text-[#9E9893]" text='BigM' onClick={() => copyToClipboard(bigM)} />
        <Card classname="bg-[#584F44] text-[#9E9893]" text='simplex' onClick={() => copyToClipboard(simplex)} />
        <Card classname="bg-[#584F44] text-[#9E9893]" text='DualSimplex' onClick={() => copyToClipboard(dualSimplex)} />
        <Card classname="bg-[#584F44] text-[#9E9893]" text='LCM' onClick={() => copyToClipboard(LCM)} />
        <Card classname="bg-[#584F44] text-[#9E9893]" text='SteepestDescent' onClick={() => copyToClipboard(SD)} />
      </div>
    </div>
  );
}

export default C2C;
