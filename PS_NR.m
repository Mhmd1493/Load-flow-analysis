% Line parameters
Z12 = 0.01008 + 0.05040i;
Y12 = 1 / Z12;
Ysh12 = 0.05125 / 2;

Z13 = 0.00744 + 0.03720i;
Y13 = 1 / Z13;
Ysh13 = 0.03875 / 2;

Z24 = 0.00744 + 0.03720i;
Y24 = 1 / Z24;
Ysh24 = 0.03875 / 2;

Z34 = 0.01272 + 0.06360i;
Y34 = 1 / Z34;
Ysh34 = 0.06375 / 2;

% Base values
Sbase = 100e6; % 100 MVA
Vbase = 230e3; % 230 kV

% Bus data
Pgen1 = 0;
Qgen1 = 0;
Pload1 = 50e6;
Qload1 = 30.99e6;
V1 = 1.0;

Pgen2 = 0;
Qgen2 = 0;
Pload2 = 170e6;
Qload2 = 105.35e6;
V2 = 1.0;

Pgen3 = 0;
Qgen3 = 0;
Pload3 = 170e6;
Qload3 = 105.35e6;
V3 = 1.0;

Pgen4 = 318e6;
Qgen4 = 0;
Pload4 = 80e6;
Qload4 = 49.58e6;
V4 = 1.02;

% Initialize variables
tolerance = 1e-6;
maxIterations = 20;

% Initialize bus voltages
V = [V1; V2; V3; V4];

% Initialize power injections
P = [Pgen1 - Pload1; Pgen2 - Pload2; Pgen3 - Pload3; Pgen4 - Pload4] / Sbase;
Q = [Qgen1 - Qload1; Qgen2 - Qload2; Qgen3 - Qload3; Qgen4 - Qload4] / Sbase;

% Initialize Jacobian matrix
J = zeros(4);

% Initialize mismatch vector
mismatch = ones(4, 1);

% Start iteration
iteration = 0;

while max(abs(mismatch)) > tolerance && iteration < maxIterations
    iteration = iteration + 1;
    
    % Calculate admittance matrix elements
    G = [real(Y12) + real(Y13), -real(Y12), -real(Y13), 0;
         -real(Y12), real(Y12) + real(Y24), 0, -real(Y24);
         -real(Y13), 0, real(Y13) + real(Y34), -real(Y34);
         0, -real(Y24), -real(Y34), real(Y24) + real(Y34)];
    
    B = [imag(Y12) + imag(Y13) + Ysh12 + Ysh13, -imag(Y12), -imag(Y13), 0;
         -imag(Y12), imag(Y12) + imag(Y24) + Ysh12 + Ysh24, 0, -imag(Y24);
         -imag(Y13), 0, imag(Y13) + imag(Y34) + Ysh13 + Ysh34, -imag(Y34);
         0, -imag(Y24), -imag(Y34), imag(Y24) + imag(Y34) + Ysh24 + Ysh34];
    
    % Calculate power injections and mismatches
    S = V .* conj(G * V) + V .* conj(B * V);
    mismatch = [P - real(S); Q - imag(S)];
    
    % Calculate Jacobian matrix
    for i = 1:4
        for j = 1:4
            if i == j
                J(i, j) = -imag(Y12 + Y13 + Ysh12 + Ysh13) * abs(V(i))^2 - Q(i) + imag(Y12) * abs(V(2))^2 ...
                    - imag(Y13) * abs(V(3))^2 + imag(Y12) * abs(V(1))^2 + imag(Y13) * abs(V(1))^2;
            else
                J(i, j) = -imag(Y12 + Y13 + Ysh12 + Ysh13) * V(i) * conj(V(j));
            end
        end
    end
    
    % Solve for voltage increments
    deltaV = J \ (-mismatch);
    
    % Update bus voltages
    V = V + deltaV;
end

% Calculate power flows
S12 = V(1) * conj(Y12 * V(1) - Y12 * V(2));
S13 = V(1) * conj(Y13 * V(1) - Y13 * V(3));
S24 = V(2) * conj(Y24 * V(2) - Y24 * V(4));
S34 = V(3) * conj(Y34 * V(3) - Y34 * V(4));

% Calculate power losses
Ploss12 = real(S12);
Ploss13 = real(S13);
Ploss24 = real(S24);
Ploss34 = real(S34);

% Convert results to per unit
V = abs(V);
P = P * Sbase;
Q = Q * Sbase;
Ploss12 = Ploss12 * Sbase;
Ploss13 = Ploss13 * Sbase;
Ploss24 = Ploss24 * Sbase;
Ploss34 = Ploss34 * Sbase;