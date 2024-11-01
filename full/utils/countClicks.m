function C = countClicks(Nl, D0, D1, l, a, b)
    %countClicks: This function counts the occurrences of different detection 
    %              click events (no clicks, clicks on D0, clicks on D1, and clicks 
    %              on both) across `Nl` pulse intensities for both matching and 
    %              non-matching basis choices.
    %
    % Inputs:
    %     Nl  - Number of pulse intensities
    %     D0  - Detection event at detector 0 (binary)
    %     D1  - Detection event at detector 1 (binary)
    %     l   - Intensity level of each pulse
    %     a   - Alice's basis choice (binary)
    %     b   - Bob's basis choice (binary)
    %
    % Outputs:
    %     C   - A 1 x (8 * Nl) array that contains the counts of each detection
    %           event for all intensities. The first half of the array corresponds 
    %           to non-matching bases (a ~= b) and the second half corresponds to 
    %           matching bases (a == b).
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    % Initialize the count matrix for Nl intensities, 4 detection events, and 2 basis matching conditions
    C = zeros(Nl, 4, 2);
    
    % Determine if the basis choices match or not
    match = (a == b) + 1;  % 1 for non-matching, 2 for matching

    % Loop through each pulse intensity
    for i = 1:Nl
        % Loop for both matching (m=2) and non-matching (m=1) cases
        for m = 1:2
            % Count the four types of click events for the current intensity and match condition
            C(i, 1, m) = sum(~D0 & ~D1 & l == i & m == match);  % No clicks (00)
            C(i, 2, m) = sum( D0 & ~D1 & l == i & m == match);  % Click on D0 (01)
            C(i, 3, m) = sum(~D0 &  D1 & l == i & m == match);  % Click on D1 (10)
            C(i, 4, m) = sum( D0 &  D1 & l == i & m == match);  % Click on both (11)
        end
    end
    
    % Reshape the count matrix into a 1 x (8 * Nl) vector:
    % First half for non-matching, second half for matching.
    C = [reshape(C(:,:,1), 1, []) reshape(C(:,:,2), 1, [])];
end