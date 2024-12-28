function [Ps, d_thetaR] = Pabxe(thetaR, thetaF, varR, varF, a, b, x, e)
    %Pabxe: Computes the detection probabilities and their derivatives for 
    %       all intensities (lambdas), w.r.t. the parameters in thetaR. 
    %       The detection probabilities depend on Alice's and Bob's basis 
    %       choices (a, b), Alice's bit choice (x), and Eve's interception 
    %       flag (e)
    %
    % Inputs:
    %   thetaR - Values of parameters listed in varR 
    %   thetaF - Values of parameters listed in varF 
    %   varR   - Cell array of random variable names
    %   varF   - Cell array of fixed variable names
    %   a      - Binary input, Alice's basis choice (0 or 1)
    %   b      - Binary input, Bob's basis choice (0 or 1)
    %   x      - Binary input, Alice's bit choice (0 or 1)
    %   e      - Binary input, Eve's interception flag 
    %            (1 if interception occurs, 0 otherwise)
    %
    % Outputs:
    %   Ps       - Combined detection probabilities for all intensities (lambdas)
    %   d_thetaR - Cell array of partial derivatives of the detection probabilities
    %              with respect to the parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Concatenate random and fixed variable lists and values
    vars = [varR(:); varF(:)];
    thetas = [thetaR(:); thetaF(:)];

    % Convert theta arrays into a structured format
    S = cell2struct(thetas, vars);

    % Extract structured variables
    lambdas = S.lambdas;
    alpha = S.alpha;
    dAB = S.dAB;
    pc0 = S.pc0;
    pc1 = S.pc1;
    pd0 = S.pd0;
    pd1 = S.pd1;
    pe = S.pe;
    dAE = S.dAE;
    pEB = S.pEB;
    k = S.k;

    % Number of intensities
    Nl = size(lambdas, 2);

    % Compute transmission probabilities based on distances and attenuation
    pAB = dist2prob(dAB, alpha);
    pAE = dist2prob(dAE, alpha);

    % % Calculate beam-splitter probabilities and their derivatives
    [p0, dp0_pe] = p_bs(0, x, a, b, pe);
    p1 = 1 - p0;
    dp1_pe = -dp0_pe;

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

    % Set up partial derivatives placeholders
    hat_pc = {hat_pc_or, hat_pc_x1, hat_pc_1x};
    hat_pd = {pd_or, pd0, pd1};
    hat_k = e * k;

    % Preallocate the output probabilities
    Ps = cell(1, 3);

    % Compute partial derivatives if requested
    if nargout > 1
        d_ = struct();
        outs = [true, false, false, false, false];  % Default: only Ps is computed

        % Check for variables in varR related to lambda
        [~, hat_lambda_inds] = ismember({'lambdas', 'alpha', 'dAB', 'dAE'}, varR);

        if any(hat_lambda_inds)
            dP_d_hat_lambda = cell(1, 3);

            % Compute lambda-related derivatives
            if hat_lambda_inds(1)
                d_.lambdas = cell(1, 3);
                d_hat_lambda_d_lambda = (pAE.^e) .* (pAB.^(1 - e));
            end
            if hat_lambda_inds(2)
                d_.alpha = cell(1, 3);
                da_dAB = -log(10) * dAB .* pAB / 10;
                da_dAE = -log(10) * dAE .* pAE / 10;
                d_hat_lambda_d_alpha = lambdas .* (da_dAE.^e) .* (da_dAB.^(1 - e));
            end
            if hat_lambda_inds(3)
                d_.dAB = cell(1, 3);
                dd_dAB = -log(10) * alpha .* pAB / 10;
                d_hat_lambda_d_dAB = lambdas .* dd_dAB * (1 - e);
            end
            if hat_lambda_inds(4)
                d_.dAE = cell(1, 3);
                dd_dAE = -log(10) * alpha .* pAE / 10;
                d_hat_lambda_d_dAE = lambdas .* dd_dAE * e;
            end
            outs(2) = true;
        end

        % Check for variables in varR related to pc
        [~, hat_pc_inds] = ismember({'pc0', 'pc1', 'pe', 'pEB'}, varR);

        if any(hat_pc_inds)
            dP_d_hat_pc = cell(1, 3);

            % Compute pc-related derivatives
            if hat_pc_inds(1)
                d_.pc0 = cell(1, 3);
                d_hat_pc_d_pc0 = {p0 .* pEB.^e, p0 .* pEB.^e, 0};
            end
            if hat_pc_inds(2)
                d_.pc1 = cell(1, 3);
                d_hat_pc_d_pc1 = {p1 .* pEB.^e, 0, p1 .* pEB.^e};
            end
            if hat_pc_inds(3)
                d_.pe = cell(1, 3);
                d_hat_pc_d_pe = {-pEB.^e .* (pc1 - pc0) .* dp0_pe, pEB.^e .* pc0 .* dp0_pe, pEB.^e .* pc1 .* dp1_pe};
            end
            if hat_pc_inds(4)
                d_.pEB = cell(1, 3);
                d_hat_pc_d_pEB = {pc_or .* e, pc0 .* p0 .* e, pc1 .* p1 .* e};
            end
            outs(3) = true;
        end

        % Check for variables in varR related to pd
        [~, hat_pd_inds] = ismember({'pd0', 'pd1'}, varR);

        if any(hat_pd_inds)
            dP_d_hat_pd = cell(1, 3);

            % Compute pd-related derivatives
            if hat_pd_inds(1)
                d_.pd0 = cell(1, 3);
                d_hat_pd_d_pd0 = {qd1, 1, 0};
            end
            if hat_pd_inds(2)
                d_.pd1 = cell(1, 3);
                d_hat_pd_d_pd1 = {qd0, 0, 1};
            end
            outs(4) = true;
        end

        % Check for variables in varR related to k
        [~, hat_hat_ind] = ismember('k', varR);
        if hat_hat_ind
            dP_d_hat_k = cell(1, 3);
            d_.k = cell(1, 3);
            d_hat_k_d_k = e;
            outs(5) = true;
        end

        % Precompute zero derivatives for fixed variables
        d_zeros = zeros(size(lambdas));

        if ismember('pa0', varR)
            d_.pa0 = {d_zeros, d_zeros, d_zeros};
        end
        if ismember('pa1', varR)
            d_.pa1 = {d_zeros, d_zeros, d_zeros};
        end
        if ismember('Delta', varR)
            d_.Delta = {d_zeros, d_zeros, d_zeros};
        end

        % Compute probabilities and derivatives
        for i = 1:3 % for each case (or, x1, 1x)
            varout = cell(1, sum(outs));
            [varout{:}] = PNS(hat_lambda, hat_pc{i}, hat_pd{i}, hat_k, outs);

            Ps{i} = varout{1};  % Store probability

            % Compute partial derivatives w.r.t. lambdas
            if outs(2)
                dP_d_hat_lambda{i} = varout{sum(outs(1:2))};
                if hat_lambda_inds(1), d_.lambdas{i} = dP_d_hat_lambda{i} .* d_hat_lambda_d_lambda; end
                if hat_lambda_inds(2), d_.alpha{i} = dP_d_hat_lambda{i} .* d_hat_lambda_d_alpha; end
                if hat_lambda_inds(3), d_.dAB{i} = dP_d_hat_lambda{i} .* d_hat_lambda_d_dAB; end
                if hat_lambda_inds(4), d_.dAE{i} = dP_d_hat_lambda{i} .* d_hat_lambda_d_dAE; end
            end

            % Compute partial derivatives w.r.t. pc
            if outs(3)
                dP_d_hat_pc{i} = varout{sum(outs(1:3))};
                if hat_pc_inds(1), d_.pc0{i} = dP_d_hat_pc{i} .* d_hat_pc_d_pc0{i}; end
                if hat_pc_inds(2), d_.pc1{i} = dP_d_hat_pc{i} .* d_hat_pc_d_pc1{i}; end
                if hat_pc_inds(3), d_.pe{i} = dP_d_hat_pc{i} .* d_hat_pc_d_pe{i}; end
                if hat_pc_inds(4), d_.pEB{i} = dP_d_hat_pc{i} .* d_hat_pc_d_pEB{i}; end
            end

            % Compute partial derivatives w.r.t. pd
            if outs(4)
                dP_d_hat_pd{i} = varout{sum(outs(1:4))};
                if hat_pd_inds(1), d_.pd0{i} = dP_d_hat_pd{i} .* d_hat_pd_d_pd0{i}; end
                if hat_pd_inds(2), d_.pd1{i} = dP_d_hat_pd{i} .* d_hat_pd_d_pd1{i}; end
            end

            % Compute partial derivatives w.r.t. k
            if outs(5)
                dP_d_hat_k{i} = varout{sum(outs)};
                d_.k{i} = dP_d_hat_k{i} .* d_hat_k_d_k;
            end
        end

        % Combine partial derivatives into d_thetaR
        d_thetaR = cell(1, numel(varR));
        for i = 1:length(d_thetaR)
            d_thetaR{i} = d_.(varR{i});
        end

        % Derive d_thetaR w.r.t. joint probabilities
        for i = 1:numel(d_thetaR)
            d_thetaR00 = -d_thetaR{i}{1};
            d_thetaR01 =  d_thetaR{i}{1} - d_thetaR{i}{3};
            d_thetaR10 =  d_thetaR{i}{1} - d_thetaR{i}{2};
            d_thetaR11 =  d_thetaR{i}{2} + d_thetaR{i}{3} - d_thetaR{i}{1};
            d_thetaR{i} = [d_thetaR00, d_thetaR01, d_thetaR10, d_thetaR11] / Nl;
        end

        % Special case: handle lambdas separately
        if hat_lambda_inds(1)
            d_lambdas = cell(1, Nl);
            for i = 1:Nl
                idx = zeros(1, Nl * 4);
                idx((i - 1) + 1:Nl:Nl * 4) = 1;
                d_lambdas{i} = d_thetaR{hat_lambda_inds(1)} .* idx;
            end
            d_thetaR = [d_thetaR(1:(hat_lambda_inds(1) - 1)), d_lambdas(:)', d_thetaR((hat_lambda_inds(1) + 1):end)];
        end

    else
        % If no derivatives are requested, only compute detection probabilities
        for i = 1:3 % for each case (or, x1, 1x)
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