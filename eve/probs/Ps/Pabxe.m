function [Ps, d_thetaE] = Pabxe(thetaA, thetaB, thetaE, a, b, x, e)
    %Pabxe: This function computes the detection probabilities at Bob's detectors 
    %       for given system parameters. It returns detection probabilities for 
    %       four possible outcomes at Bob’s two detectors (00, 01, 10, 11) and 
    %       calculates partial derivatives with respect to parameters in thetaE.
    %
    % Inputs:
    %   thetaA - Cell array of system parameters specific to the channel between 
    %            Alice and Bob, containing [lambdas, alpha, dAB].
    %   thetaB - Cell array of system parameters for Bob's detectors, containing 
    %            [pa0, pa1, pc0, pc1, pd0, pd1, pe].
    %   thetaE - Cell array of system parameters related to Eve’s eavesdropping, 
    %            containing [dAE, pEB, k].
    %   a      - Binary input, Alice's basis choice for the pulse (0 or 1).
    %   b      - Binary input, Bob's basis choice for the pulse (0 or 1).
    %   x      - Binary input, Alice's bit choice for the pulse (0 or 1).
    %   e      - Binary input, Eve's interception flag 
    %            (e=1 if she intercepted the pulse, e=0 otherwise).
    %
    % Outputs:
    %   Ps       - Array of detection probabilities at Bob’s two detectors for the 
    %              four possible detection outcomes [P_00, P_01, P_10, P_11], 
    %              normalized by the number of lambda values (Nl).
    %   d_thetaE - Cell array of partial derivatives of detection probabilities 
    %              with respect to the parameters in thetaE.
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Extract parameters from input arrays
    [lambdas, alpha, dAB] = deal(thetaA{:});
    [~, ~, pc0, pc1, pd0, pd1, pe] = deal(thetaB{:});
    [dAE, pEB, k, ~] = deal(thetaE{:});

    % Extract number of intensity values
    Nl = length(lambdas);

    % Compute transmission probabilities based on distances between Alice, Bob, and Eve
    pAB = dist2prob(dAB, alpha);
    pAE = dist2prob(dAE, alpha);

    % Calculate beam-splitter probabilities and complementary probabilities
    p0 = p_bs(0, x, a, b, pe);
    p1 = 1 - p0;

    qd0 = 1 - pd0;
    qd1 = 1 - pd1;

    % Adjust lambda based on interception probabilities
    hat_lambda = lambdas .* pAB.^(1 - e) .* pAE.^e;

    % Calculate mixed and effective detection probabilities at Bob's detectors
    pc_mix = pc1 - p0 .* (pc1 - pc0);
    hat_pc_or = pc_mix .* pEB.^e;
    hat_pc_x1 = p0 .* pc0 .* pEB.^e;
    hat_pc_1x = p1 .* pc1 .* pEB.^e;
    hat_pd_or = 1 - qd0 .* qd1;

    % Aggregate detection probabilities into cells for processing
    hat_pc = {hat_pc_or, hat_pc_x1, hat_pc_1x};
    hat_pd = {hat_pd_or, pd0, pd1};
    hat_k = e * k;
    Ps = cell(1, 3);

    % If partial derivatives are requested
    if nargout > 1
        % Define requested outputs (ignoring only d_pd)
        outs = [true, true, true, false, true];
        
        % Derivative of hat_lambda w.r.t. dAE
        d_dAE = cell(1, 3);
        dd_dAE = -log(10) * alpha .* pAE / 10;
        d_hat_lambda_d_dAE = lambdas .* dd_dAE * e;

        % Derivative of hat_pc w.r.t. pEB
        d_pEB = cell(1, 3);
        d_hat_pc_d_pEB = {pc_mix .* e, pc0 .* p0 .* e, pc1 .* p1 .* e};

        % Derivative of hat_k w.r.t. k
        d_k = cell(1, 3);
        d_hat_k_d_k = e;

        % Calculate probabilities and derivatives for each case (1: or, 2: x1, and 3: 1x)
        for i = 1:3
            [Ps{i}, dP_d_hat_lambda, dP_d_hat_pc, dP_d_hat_k] = PNS(hat_lambda, hat_pc{i}, hat_pd{i}, hat_k, outs);

            % Calculate derivatives of Por, Px1, or P1x with respect to parameters in
            % thetaE (only dAE, pEB, and k are involved)
            d_dAE{i} = dP_d_hat_lambda .* d_hat_lambda_d_dAE;
            d_pEB{i} = dP_d_hat_pc .* d_hat_pc_d_pEB{i};
            d_k{i} = dP_d_hat_k .* d_hat_k_d_k;
        end

        % Precompute zero derivatives for fixed variables
        d_zeros = zeros(size(lambdas));
        d_Delta = {d_zeros, d_zeros, d_zeros};

        % Derive the derivatives of P00, P01, P10, and P11 from the
        % derivatives of Por, Px1, and P1x and store them into d_thetaE
        d_thetaE = {d_dAE, d_pEB, d_k, d_Delta};
        for i = 1:numel(d_thetaE)
            d_thetaE00 = -d_thetaE{i}{1};
            d_thetaE01 = d_thetaE{i}{1} - d_thetaE{i}{3};
            d_thetaE10 = d_thetaE{i}{1} - d_thetaE{i}{2};
            d_thetaE11 = d_thetaE{i}{2} + d_thetaE{i}{3} - d_thetaE{i}{1};
            d_thetaE{i} = [d_thetaE00, d_thetaE01, d_thetaE10, d_thetaE11] / Nl;
        end

    else
        % Compute detection probabilities if no derivatives are requested
        for i = 1:3
            Ps{i} = PNS(hat_lambda, hat_pc{i}, hat_pd{i}, hat_k);
        end
    end

    P_or = Ps{1};
    P_x1 = Ps{2};
    P_1x = Ps{3};

    % Derive P00, P01, P10, and P11 from P_or, P_x1, and P_1x
    P_00 = 1 - P_or;
    P_01 = P_or - P_1x;
    P_10 = P_or - P_x1;
    P_11 = P_x1 + P_1x - P_or;

    % Normalize probabilities by the number of lambda values
    Ps = [P_00, P_01, P_10, P_11] / Nl;
end