%Prompt user for number of Joints and Members
J = input('Enter number of joints: ');
M = input('Enter number of members: ');

%Initialize Matrices
X = zeros(1,J);
Y = zeros(1,J);
C = zeros(J,M);
Sx = zeros(J,3);
Sy = zeros(J,3);
L = zeros(2*J,1);
A = zeros(2*J,M+3);
lengths = zeros(1,M);

%Prompt user for Joint Coordinates
for i = 1:J
    coordsPrompt = sprintf('Enter (x,y) coordinates of Joint %d in the form "x y": ', i);
    coords = "" + split(input(coordsPrompt, 's'));
    xcoord = str2num(coords(1));
    ycoord = str2num(coords(2));
    X(i) = xcoord;
    Y(i) = ycoord;
end

%Prompt user for which Joints connect which Members  
joint1 = 0;
joint2 = 0;
for i = 1:M
    jointPrompt = sprintf('Enter number of the Joints connected to Member %d in the form "# #": ', i);
    joints = "" + split(input(jointPrompt, 's'));
    joint1 = str2num(joints(1));
    joint2 = str2num(joints(2));
    C(joint1,i) = 1;
    C(joint2,i) = 1;
    lengths(i) = sqrt((X(joint2) - X(joint1))^2 + (Y(joint2) - Y(joint1))^2);
end

%Prompt user for reaction forces
numRxnX = input('Enter number of reaction forces in the x-dir: \n');
count = 0;
for i = 1:numRxnX
    jointNumPrompt = sprintf('Enter number of Joint reaction force %d acts on: ', i);
    jointNum = input(jointNumPrompt);
    Sx(jointNum,1) = 1;
    count = count + 1;
end
numRxnY = input('Enter number of reaction forces in the y-dir: \n');
for i = 1:numRxnY
    jointNumPrompt = sprintf('Enter number of Joint reaction force %d acts on: ', i);
    jointNum = input(jointNumPrompt);
    Sy(jointNum,count + i) = 1;
end

%Save Matrices to File
save('TrussDesign1_CharlesDominicBharath_A1.mat','C','Sx','Sy','X','Y','L');

%Set up A matrix
for r1 = 1:J
    firstx = 0;
    secondx = 0;
    firsty = 0;
    secondy = 1;
    for c = 1:M
        firstx = 0;
        firsty = 0;
        secondx = 0;
        secondy = 1;
        if (C(r1,c) == 1)
            firstx = X(r1);
            firsty = Y(r1);
            for r2 = 1:J
                if (C(r2,c) == 1)
                    if (r2 ~= r1)
                        secondx = X(r2);
                        secondy = Y(r2);
                    end
                end
            end
        end
        xUnitVector = secondx - firstx;
        distance = sqrt( (secondx - firstx)^2 +  (secondy - firsty)^2 );
        A(r1,c) = xUnitVector/distance;
    end    
end

for r1 = 1:J
    firstx = 0;
    secondx = 1;
    firsty = 0;
    secondy = 0;
    for c = 1:M
        firstx = 0;
        firsty = 0;
        secondx = 1;
        secondy = 0;
        if (C(r1,c) == 1)
            firstx = X(r1);
            firsty = Y(r1);
            for r2 = 1:J
                if (C(r2,c) == 1)
                    if (r2 ~= r1)
                        secondx = X(r2);
                        secondy = Y(r2);
                    end
                end
            end
        end
        yUnitVector = secondy - firsty;
        distance = sqrt( (secondx - firstx)^2 +  (secondy - firsty)^2 );
        A(r1+J,c) = yUnitVector/distance;
    end  
end
[r, c] = size(A);
A(1:J,M+1:c) = Sx;
A(J+1:r,M+1:c) = Sy;

%Set up L matrix
loadJointNum = input('Enter joint # that load is applied to: \n');
loadWeight = input('Enter load weight in oz: \n');
L(J + loadJointNum) = loadWeight;

%Solve for T
T=linsolve(A,L);

%Print Analysis
fprintf('EK301, Section A1, Charles P., Domininc J., Bharath V., 4/1/2023\n')
fprintf('Load: %d oz\n', loadWeight);
fprintf('Member Forces in oz: \n');
for i = 1:M
    comptens = "";
    if (T(i) < 0)
        comptens = comptens + "(Compression)";
    else (T(i) > 0);
        comptens = comptens + "(Tension)";
    end
    fprintf('m%d: %.3f %s\n', i, abs(T(i)), comptens);
end
fprintf('Reaction Forces in oz:\n')
for i = 1:numRxnX
    fprintf('Sx%d: %.3f\n', i, T(i+M));
end
for i = 1:numRxnY
    fprintf('Sy%d: %.3f\n', i, T(i+numRxnX+M));
end
cost = 10*J + sum(lengths);
fprintf('Cost of truss: %.2f\n', cost);


Pc = zeros(1,length(T)-3);
R = zeros(1,length(T)-3);
for i = 1:length(T)-3
    Pc(i)=4338/(lengths(i)^2.125);
    R(i)=T(i)/Pc(i);
end
Wfailure = -Pc./R;
Memory = min(Wfailure(Wfailure>0));
bMem = zeros(1,M);
in = zeros(1,M);
for i=length(Wfailure)
    if Memory == Wfailure(1,i)
        bMem(1,i) =i;
        in(1,i) = i;
    end
end
minIndex=1;
for i=2:M
    if (Wfailure(i) < Wfailure(minIndex))
        minIndex=i;
    end
end

magicfind = find(R(1,:) == min(R));
loadmax = - loadWeight/ R(1, magicfind(1));
error = 1.01 * loadmax / Pc(magicfind(1));
fprintf('Theoretical max load/cost ratio in oz/$: %.4f\n', loadmax/cost);
fprintf('The Critical Member is Member %d\n', minIndex);
fprintf('Length of Member %d: %.3f in\n', minIndex, lengths(minIndex));
fprintf('Predicted Buckling Strength of Member %d: %.3f +/- %.3f oz\n', minIndex, Pc(minIndex), Pc(minIndex)*0.05);

fprintf('The theoretical Max Load is: %.3f +/- %.3f oz \n', loadmax, error);
