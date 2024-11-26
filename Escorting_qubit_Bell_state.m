% Escorting Qubit Method
%% Global definition of symbolic variables

clc; clear;
syms r q al be r1 r2
assume((0 <= r) & (r <= 1))
assume((0 <= r1) & (r1 <= 1))
assume((0 <= r2) & (r2 <= 1))
assume((0 <= q) & (q <= 1))
assume((0 <= al) & (al <= 1))
assume((0 <= be) & (be <= 1))
assume(r, 'real')
assume(r1, 'real')
assume(r2, 'real')
assume(q, 'real')
assume(al, 'real')
assume(be, 'real')

%% Define intial state

% Define Bell state
H = [1; 0]; % |0>
V = [0; 1]; % |1>
bell_state = 1/sqrt(2) * (kron(H, H) +  kron(V, V)); % Bell state 
disp('Bell state:');
disp(bell_state);

% Define ancilla states
ancilla1 = [0; 1]; 
ancilla2 = [0; 1]; 

% Initial state: Bell state  anc1  anc2
initial_state = kron(bell_state, kron(ancilla1, ancilla2));


% Density matrix of the initial state
rho_bell0 = initial_state * initial_state';


%% The First CNOT gates

% Define CNOT gates
C13 = [1 0 0 0 0 0 0 0;
   0 1 0 0 0 0 0 0;
   0 0 1 0 0 0 0 0;
   0 0 0 1 0 0 0 0;
   0 0 0 0 0 1 0 0;
   0 0 0 0 1 0 0 0;
   0 0 0 0 0 0 0 1;
   0 0 0 0 0 0 1 0]; % Three-qubit CNOT: 1 as control, 3 as target
CNOT_1 = kron(C13, eye(2)); 
CNOT_2 = kron(eye(2), C13);

% Apply CNOT gates
rho_c = CNOT_2 * CNOT_1 * rho_bell0 * CNOT_1' * CNOT_2';

%% ADC

% Define ADC Kraus operators
E0_r1 = [1, 0; 0, sqrt(1 - r1)];
E1_r1 = [0, sqrt(r1); 0, 0];
K1 = {E0_r1, E1_r1};   %q_1 & ancilla1)

E0_r2 = [1, 0; 0, sqrt(1 - r2)];
E1_r2 = [0, sqrt(r2); 0, 0];
K2 = {E0_r2, E1_r2};   %q_2 & ancilla2)

% Apply ADC to the state
rho_bell_adc = zeros(size(rho_c));
for i1 = 1:2
    for i2 = 1:2
        for j1 = 1:2
            for j2 = 1:2
                K = kron(K1{i1}, kron(K2{j1}, kron(K1{i2}, K2{j2})));
                rho_bell_adc = rho_bell_adc + K * rho_c * K';
            end
        end
    end
end

disp('Density matrix after ADC:');
disp(simplify(rho_bell_adc));

%% The Second CNOT gates

% Apply second CNOT gates
rho_final = CNOT_2 * CNOT_1 * rho_bell_adc * CNOT_1' * CNOT_2';
disp('Density matrix after second CNOT gates:');
disp(simplify(rho_final));

%% Measurement

% Define measurement operator for |11>
M11 = kron(eye(4), kron([0 0; 0 1], [0 0; 0 1])); 

% Measure |11
rho_anc_11 = simplify(M11 * rho_final * M11');
p_3 = trace(rho_anc_11); % Probability of |11>


% Display results
disp('|11> probability:');
disp(p_3);
disp('|11> density matrix:');
disp(rho_anc_11);
% rho_anc_11 = simplify(rho_anc_11 / p_3); % Normalize density matrix
% disp('|11 density matrix (Normalized):');
% disp(rho_anc_11);
%% Final

% Partial trace over ancilla qubits
disp('Partial trace over ancillas:');
reduced_rho = PartialTrace(rho_anc_11, [3, 4]); % Trace out ancilla qubits
disp(simplify(reduced_rho));

disp('Output normalized state:');
rout=reduced_rho/p_3;
disp(simplify(rout));


%%
function [result] = PartialTrace(rho, subsys)
% PartialTrace performs a partial trace operation on the density matrix rho
% over the specified subsystem indices subsys.
% 
% Parameters:
%   rho: Joint system density matrix (2^N x 2^N)
%   subsys: Indices of subsystems to trace out (e.g., [3,4])
% 
% Returns:
%   result: Density matrix after tracing out the specified subsystems

% Determine the total number of qubits in the joint system
    total_qubits = log2(size(rho, 1));

% Ensure the provided subsystem indices are valid
    if any(subsys > total_qubits) || any(subsys < 1)
        error('Subsystem indices exceed the range of the joint system');
    end

% Determine the Hilbert space dimension of the subsystems to trace out
    trace_qubits = length(subsys);
    trace_dim = 2^trace_qubits;

% Determine the subsystems to retain
    keep_subsys = setdiff(1:total_qubits, subsys);
    keep_dim = 2^length(keep_subsys);

% Initialize the resulting density matrix
    result = zeros(keep_dim, keep_dim);
    
% Iterate over the standard basis vectors of the traced-out subsystems
    for i = 0:(trace_dim - 1)
% Generate the binary representation of the basis vector index
        idx = dec2bin(i, trace_qubits) - '0';

% Create the basis state for the traced-out subsystems
        state = 1; % Initialize
            for k = 1:trace_qubits
                if idx(k) == 0
                   state = kron(state, [1; 0]); % Basis state |0
                else
                   state = kron(state, [0; 1]); % Basis state |1
                end
            end

% Expand the basis state to the full joint system
        full_state = 1; % Initialize as identity
        subsys_idx = 1;
    for q = 1:total_qubits
        if ismember(q, subsys) % Subsystems to trace out
            if idx(subsys_idx) == 0
                full_state = kron(full_state, [1; 0]); % |0
            else
                full_state = kron(full_state, [0; 1]); % |1
            end
                subsys_idx = subsys_idx + 1;
            else
                full_state = kron(full_state, eye(2)); % Retained subsystems use identity
            end
    end

% Perform the partial trace for the current basis vector
        result = result + full_state' * rho * full_state;
    end
 end