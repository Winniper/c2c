import Card from '../components/Card'

const C2C = () => {
  const graphCode = `\
% Graphical Method for Linear Programming
format short
clear all
clc

% Define constraints
a = [1 2; 1 1; 1 -2];
b = [10; 6; 1];
c = [5 3];

x1 = 0 : max(b);
x21 = (b(1) - a(1,1)*x1) / a(1,2);
x22 = (b(2) - a(2,1)*x1) / a(2,2);
x23 = (b(3) - a(3,1)*x1) / a(3,2);

x21 = max(0, x21);
x22 = max(0, x22);
x23 = max(0, x23);

% Plot constraints
plot(x1, x21, 'r', x1, x22, 'k', x1, x23, 'p');
hold on;

% Finding corner points
cx1 = find(x1 == 0);
c1 = find(x21 == 0);
line1 = [x1([c1, cx1]); x21([c1, cx1])]';

c2 = find(x22 == 0);
line2 = [x1([c2, cx1]); x22([c2, cx1])]';

c3 = find(x23 == 0);
line3 = [x1([c3, cx1]); x23([c3, cx1])]';

cornpt = unique([line1; line2; line3], 'rows');

% Intersecting points of constraints
pt = [0; 0];
for i = 1 : size(a,1)
    for j = i+1 : size(a,1)
        A1 = a([i j], :);
        B1 = b([i j], :);
        X = inv(A1) * B1;
        pt = [pt X];
    end
end
ptt = pt';

% All feasible corner points
allpts = [ptt; cornpt];
points = unique(allpts, 'rows');

% Check feasibility
for i = 1 : size(points,1)
    const1(i) = a(1,1)*points(i,1) + a(1,2)*points(i,2) - b(1);
    const2(i) = a(2,1)*points(i,1) + a(2,2)*points(i,2) - b(2);
    const3(i) = a(3,1)*points(i,1) + a(3,2)*points(i,2) - b(3);
end
s1 = find(const1 > 0);
s2 = find(const2 > 0);
s3 = find(const3 > 0);
S = unique([s1 s2 s3]);
points(S, :) = [];

% Calculate objective function value
value = points * c';
table = [points value];
[obj, index] = max(value);
X1 = points(index,1);
X2 = points(index,2);
fprintf('Objective value is %f at (%f,%f)', obj, X1, X2);
`;

const graphCodeWith2Constraints =`\
% Graphical Method for Linear Programming
format short
clear all
clc

% Define constraints
a = [1 2; 1 1];
b = [10; 6];
c = [5 3];

x1 = 0 : max(b);
x21 = (b(1) - a(1,1)*x1) / a(1,2);
x22 = (b(2) - a(2,1)*x1) / a(2,2);


x21 = max(0, x21);
x22 = max(0, x22);

% Plot constraints
plot(x1, x21, 'r', x1, x22, 'k');
hold on;

% Finding corner points
cx1 = find(x1 == 0);
c1 = find(x21 == 0);
line1 = [x1([c1, cx1]); x21([c1, cx1])]';

c2 = find(x22 == 0);
line2 = [x1([c2, cx1]); x22([c2, cx1])]';



cornpt = unique([line1; line2], 'rows');

% Intersecting points of constraints
pt = [0; 0];
for i = 1 : size(a,1)
    for j = i+1 : size(a,1)
        A1 = a([i j], :);
        B1 = b([i j], :);
        X = inv(A1) * B1;
        pt = [pt X];
    end
end
ptt = pt';

% All feasible corner points
allpts = [ptt; cornpt];
points = unique(allpts, 'rows');

% Check feasibility
for i = 1 : size(points,1)
    const1(i) = a(1,1)*points(i,1) + a(1,2)*points(i,2) - b(1);
    const2(i) = a(2,1)*points(i,1) + a(2,2)*points(i,2) - b(2);

end
s1 = find(const1 > 0);
s2 = find(const2 > 0);

S = unique([s1 s2]);
points(S, :) = [];

% Calculate objective function value
value = points * c';
table = [points value];
[obj, index] = max(value);
X1 = points(index,1);
X2 = points(index,2);
fprintf('Objective value is %f at (%f,%f)', obj, X1, X2);
`;

  const bfsCode = `\
format short
clear all
clc

% Basic Feasible Solution Method for Linear Programming

A = [2 3 -1 4; 1 2 6 -7];
b = [8; 3];
c = [2 3 4 7];

% Number of equations and variables
m = size(A,1);
n = size(A,2);

if n >= m
    nCm = nchoosek(n, m);
    P = nchoosek(1:n, m);
    sol = [];

    for i = 1 : nCm
        y = zeros(n,1);
        A1 = A(:, P(i,:));
        X = inv(A1) * b;

        if all(X >= 0 & X ~= inf & X ~= -inf)
            y(P(i,:)) = X;
            sol = [sol y];
        end
    end

    % Compute objective function values
    z = c * sol;

    % Display the results
    fprintf("Basic Feasible Solutions and Objective Function Values:\\n");
    fprintf("------------------------------------------------------\\n");
    fprintf("  x1      x2      x3      x4      Z (Objective)\\n");
    fprintf("------------------------------------------------------\\n");

    for i = 1:size(sol,2)
        fprintf("%6.2f  %6.2f  %6.2f  %6.2f  %10.2f\\n", sol(:,i)', z(i));
    end
else
    disp("Not enough variables for BFS.");
end
`;

  const copyToClipboard = (text: string) => {
    navigator.clipboard.writeText(text.replace(/\n/g, "\r\n"))
  };

  return (
    <div className='bg-[#0D1117] h-dvh w-auto flex justify-center items-center'>
      <div className='flex gap-5 pt-24 pb-10'>
        <Card classname="bg-[#584F44] text-[#9E9893]" text='Long' onClick={() => copyToClipboard(graphCode)} />
        <Card classname="bg-[#584F44] text-[#9E9893]" text='Long 2.0' onClick={() => copyToClipboard(graphCodeWith2Constraints)} />
        <Card classname="bg-[#584F44] text-[#9E9893]" text='Short' onClick={() => copyToClipboard(bfsCode)} />
      </div>
    </div>
  );
}

export default C2C;
