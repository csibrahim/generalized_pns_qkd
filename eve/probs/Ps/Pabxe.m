function [Ps, d_thetaE] = Pabxe(thetaA, thetaB, thetaE, a, b, x, e)
    %Pabxe: Computes the detection probabilities and their derivatives for 
    %       all intensities (lambdas), w.r.t. the parameters in thetaE. 
    %       The detection probabilities depend on Alice's and Bob's basis 
    %       choices (a, b), Alice's bit choice (x), and Eve's interception 
    %       flag (e)
    %
    % Inputs:
    %   thetaA - Alice’s parameter set [lambdas, alpha, dAB]
    %   thetaB - Bob’s parameter set [pa0, pa1, pc0, pc1, pd0, pd1, pe]
    %   thetaE - Eve’s parameter set [dAE, pEB, k, Delta]
    %   a      - Binary input, Alice's basis choice (0 or 1)
    %   b      - Binary input, Bob's basis choice (0 or 1)
    %   x      - Binary input, Alice's bit choice (0 or 1)
    %   e      - Binary input, Eve's interception flag 
    %            (1 if interception occurs, 0 otherwise)
    %
    % Outputs:
    %   Ps       - Combined detection probabilities for all intensities (lambdas)
    %   d_thetaE - Cell array of partial derivatives of the detection probabilities
    %              with respect to the parameters in thetaE
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Extract parameters from input arrays
    [lambdas, alpha, dAB] = deal(thetaA{:});
    [~, ~, pc0, pc1, pd0, pd1, pe] = deal(thetaB{:});
    [dAE, pEB, k, ~] = deal(thetaE{:});

    % Number of intensities
    Nl = size(lambdas, 2);

    % Compute transmission probabilities based on distances and attenuation
    pAB = dist2prob(dAB, alpha);
    pAE = dist2prob(dAE, alpha);

    % Calculate beam-splitter probabilities
    p0 = p_bs(0, x, a, b, pe);
    p1 = 1 - p0;

    % Complementary probabilities for dark counts
    qd0 = 1 - pd0;
    qd1 = 1 - pd1;

    % Adjust lambda based on interception probabilities
    hat_lambda = lambdas .* pAB.^(1 - e) .* pAE.^e;

    % Union probabilities for pc and pd
    pc_or = p1 .* pc1 + p0 .* pc0;
    pd_or = 1 - qd0 .* qd1;

    % Adjusted pc based on eavesdropping
    hat_pc_or = pc_or .* pEB.^e;
    hat_pc_x1 = p0 .* pc0 .* pEB.^e;
    hat_pc_1x = p1 .* pc1 .* pEB.^e;
    

    % Aggregate detection probabilities into arrays
    hat_pc = {hat_pc_or, hat_pc_x1, hat_pc_1x};
    hat_pd = {pd_or, pd0, pd1};
    hat_k = e * k;

    % Preallocate the output probabilities for each case (or, x1, 1x)
    Ps = cell(1, 3);

    % Compute partial derivatives if requested
    if nargout > 1
        % Define outputs for PNS (derivative w.r.t. pd is irrelevant)
        outs = [true, true, true, false, true];

        % Derivatives of hat_lambda w.r.t. dAE
        d_dAE = cell(1, 3);
        dd_dAE = -log(10) * alpha .* pAE / 10;
        d_hat_lambda_d_dAE = lambdas .* dd_dAE * e;

        % Derivatives of hat_pc w.r.t. pEB
        d_pEB = cell(1, 3);
        d_hat_pc_d_pEB = {pc_or .* e, pc0 .* p0 .* e, pc1 .* p1 .* e};

        % Derivatives of hat_k w.r.t. k
        d_k = cell(1, 3);
        d_hat_k_d_k = e;

        % Compute probabilities and derivatives
        for i = 1:3 % for each case (or, x1, 1x)
            [Ps{i}, dP_d_hat_lambda, dP_d_hat_pc, dP_d_hat_k] = PNS(hat_lambda, hat_pc{i}, hat_pd{i}, hat_k, outs);

            % Compute derivatives w.r.t. dAE, pEB, and k
            d_dAE{i} = dP_d_hat_lambda .* d_hat_lambda_d_dAE;
            d_pEB{i} = dP_d_hat_pc .* d_hat_pc_d_pEB{i};
            d_k{i} = dP_d_hat_k .* d_hat_k_d_k;
        end

        % Derive d_thetaE w.r.t. joint probabilities
        d_thetaE = {d_dAE, d_pEB, d_k};
        for i = 1:numel(d_thetaE)
            d_thetaE00 = -d_thetaE{i}{1};
            d_thetaE01 =  d_thetaE{i}{1} - d_thetaE{i}{3};
            d_thetaE10 =  d_thetaE{i}{1} - d_thetaE{i}{2};
            d_thetaE11 =  d_thetaE{i}{2} + d_thetaE{i}{3} - d_thetaE{i}{1};
            d_thetaE{i} = [d_thetaE00, d_thetaE01, d_thetaE10, d_thetaE11] / Nl;
        end
    else
        % If no derivatives are requested, only compute detection probabilities
        for i = 1:3  % for each case (or, x1, 1x)
            Ps{i} = PNS(hat_lambda, hat_pc{i}, hat_pd{i}, hat_k);
        end
    end

    % Extract union probabilities
    P_or = Ps{1};
    P_x1 = Ps{2};
    P_1x = Ps{3};

    % Derive joint probabilities
    P_00 = 1 - P_or;
    P_01 = P_or - P_1x;
    P_10 = P_or - P_x1;
    P_11 = P_x1 + P_1x - P_or;

    % Normalize probabilities by the number of intensities
    Ps = [P_00, P_01, P_10, P_11] / Nl;
end
