function [Ts, d_p00, d_p01, d_p10, d_p11] = Tabxel(p, q)
    %Tabxel: Computes the transition matrices for Bob's detectors and its 
    %        derivatives w.r.t. the parameters in thetaE. The detection 
    %        probabilities depend on Alice's and Bob's basis choices (a, b), 
    %        Alice's bit choice (x), Eve' interception flag (e), and 
    %        intensity (lambda).
    %
    % Inputs:
    %   p - Vector of independent click probabilities [p00, p01, p10, p11]
    %       - p00: Probability of no clicks
    %       - p01: Probability of a click on D0 only
    %       - p10: Probability of a click on D1 only
    %       - p11: Probability of clicks on both D0 and D1
    %   q - Vector of after-pulse probabilities [q0, q1]
    %       - q0: After-pulse probability for detector D0
    %       - q1: After-pulse probability for detector D1
    %
    % Outputs:
    %   Ts     - Transition matrix for the HMM
    %   d_p00  - Derivative of Ts w.r.t. p00
    %   d_p01  - Derivative of Ts w.r.t. p01
    %   d_p10  - Derivative of Ts w.r.t. p10
    %   d_p11  - Derivative of Ts w.r.t. p11
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Extract independent after-pulse probabilities
    q0 = q(1);  % After-pulse probability for detector D0
    q1 = q(2);  % After-pulse probability for detector D1
    
    % Extract independent click probabilities
    p00 = p(1);  % Probability of no clicks
    p01 = p(2);  % Probability of click on D0 only
    p10 = p(3);  % Probability of click on D1 only
    p11 = p(4);  % Probability of clicks on both D0 and D1

    % Compute marginal click probabilities
    px0 = p00 + p10;  % Probability of no click on D0
    px1 = p01 + p11;  % Probability of click on D0
    p0x = p00 + p01;  % Probability of no click on D1
    p1x = p10 + p11;  % Probability of click on D1

    % Compute joint after-pulse probabilities
    q00 = 1 - q1 - q0 + q0 * q1;  % No after-pulse
    q01 = q0 - q0 * q1;           % After-pulse on D0 only
    q10 = q1 - q0 * q1;           % After-pulse on D1 only
    q11 = q0 * q1;                % After-pulse on both detectors

    % Define the sparse transition matrix
    Ts = sparse([
                  [p00 p01 p10 p11] [0       0       0       0       0  ]
       (1 - q0) * [p00 p01 p10 p11] [p0x     0       p1x     0       0  ] * q0
       (1 - q1) * [p00 p01 p10 p11] [0       px0     0       px1     0  ] * q1
            q00 * [p00 p01 p10 p11] [p0x*q01 px0*q10 p1x*q01 px1*q10 q11]
                  [p00 p01 p10 p11] [0       0       0       0       0  ]
                  [p00 p01 p10 p11] [0       0       0       0       0  ]
       (1 - q1) * [p00 p01 p10 p11] [0       px0     0       px1     0  ] * q1
       (1 - q0) * [p00 p01 p10 p11] [p0x     0       p1x     0       0  ] * q0
                  [p00 p01 p10 p11] [0       0       0       0       0  ]
        ]);

    % If partial derivatives are requested
    if nargout > 1

        % Partial derivative w.r.t. p00
        d_p00 = sparse([
                      [1   0   0   0] [0   0   0   0   0]
           (1 - q0) * [1   0   0   0] [1   0   0   0   0] * q0
           (1 - q1) * [1   0   0   0] [0   1   0   0   0] * q1
                q00 * [1   0   0   0] [q01 q10 0   0   0]
                      [1   0   0   0] [0   0   0   0   0]
                      [1   0   0   0] [0   0   0   0   0]
           (1 - q1) * [1   0   0   0] [0   1   0   0   0] * q1
           (1 - q0) * [1   0   0   0] [1   0   0   0   0] * q0
                      [1   0   0   0] [0   0   0   0   0]
               ]);
        
        % Partial derivative w.r.t. p01
        d_p01 = sparse([
                      [0   1   0   0] [0   0   0   0   0]
           (1 - q0) * [0   1   0   0] [1   0   0   0   0] * q0
           (1 - q1) * [0   1   0   0] [0   0   0   1   0] * q1
                q00 * [0   1   0   0] [q01 0   0   q10 0]
                      [0   1   0   0] [0   0   0   0   0]
                      [0   1   0   0] [0   0   0   0   0]
           (1 - q1) * [0   1   0   0] [0   0   0   1   0] * q1
           (1 - q0) * [0   1   0   0] [1   0   0   0   0] * q0
                      [0   1   0   0] [0   0   0   0   0]
               ]);
        

        % Partial derivative w.r.t. p10
        d_p10 = sparse([
                      [0   0   1   0] [0   0   0   0   0]
           (1 - q0) * [0   0   1   0] [0   0   1   0   0] * q0
           (1 - q1) * [0   0   1   0] [0   1   0   0   0] * q1
                q00 * [0   0   1   0] [0   q10 q01 0   0]
                      [0   0   1   0] [0   0   0   0   0]
                      [0   0   1   0] [0   0   0   0   0]
           (1 - q1) * [0   0   1   0] [0   1   0   0   0] * q1
           (1 - q0) * [0   0   1   0] [0   0   1   0   0] * q0
                      [0   0   1   0] [0   0   0   0   0]
               ]);
        

        % Partial derivative w.r.t. p11
        d_p11 = sparse([
                      [0   0   0   1] [0   0   0   0   0]
           (1 - q0) * [0   0   0   1] [0   0   1   0   0] * q0
           (1 - q1) * [0   0   0   1] [0   0   0   1   0] * q1
                q00 * [0   0   0   1] [0   0   q01 q10 0]
                      [0   0   0   1] [0   0   0   0   0]
                      [0   0   0   1] [0   0   0   0   0]
           (1 - q1) * [0   0   0   1] [0   0   0   1   0] * q1
           (1 - q0) * [0   0   0   1] [0   0   1   0   0] * q0
                      [0   0   0   1] [0   0   0   0   0]
               ]);
    end
end
