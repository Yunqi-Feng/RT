function CADop_diffr = generate_edges(CADop)
r = size(CADop,1);
CADop_sca = [];
for i = 1:r
    if CADop(i,14) == 7
        CADop_sca = [CADop_sca;CADop(i,:)];
    end
end

numTriangles = size(CADop_sca, 1);

% Preallocate the output matrix: 3 edges per triangle, each edge has 6 values (2 points * 3 coords)
CADop_diffr  = zeros(numTriangles * 3, 6);

for i = 1:numTriangles
    % Extract the 3 vertices of the i-th triangle
    v1 = CADop_sca(i, 1:3);
    v2 = CADop_sca(i, 4:6);
    v3 = CADop_sca(i, 7:9);
    
    % Edge 1: v1 to v2
    CADop_diffr(3*(i-1) + 1, :) = [v1, v2];
    
    % Edge 2: v2 to v3
    CADop_diffr(3*(i-1) + 2, :) = [v2, v3];
    
    % Edge 3: v3 to v1
    CADop_diffr(3*(i-1) + 3, :) = [v3, v1];
end

zeroCoordMask = (CADop_diffr(:,1) == 0 & CADop_diffr(:,4) == 0) | ...
                (CADop_diffr(:,2) == 0 & CADop_diffr(:,5) == 0) | ...
                (CADop_diffr(:,3) == 0 & CADop_diffr(:,6) == 0);

% Keep only lines without zero coordinates
CADop_diffr = CADop_diffr(~zeroCoordMask, :);


numLines = size(CADop_diffr,1);

% Create a logical mask to mark duplicates for removal
removeMask = false(numLines,1);

% We'll use a simple nested loop for clarity (can be optimized for large datasets)
for i = 1:numLines-1
    if removeMask(i)
        continue; % already marked for removal
    end
    for j = i+1:numLines
        if removeMask(j)
            continue; % already marked for removal
        end

        CADop_diffr_A = CADop_diffr(i,:);
        CADop_diffr_B = CADop_diffr(j,:);

        % Check if lineA == lineB
        if all(CADop_diffr_A == CADop_diffr_B)
            removeMask(i) = true;
            removeMask(j) = true;

        % Check if lineA == reversed lineB
        elseif all(CADop_diffr_A(1:3) == CADop_diffr_B(4:6)) && all(CADop_diffr_A(4:6) == CADop_diffr_B(1:3))
            removeMask(i) = true;
            removeMask(j) = true;
        end
    end
end

% Remove all duplicate lines (and their duplicates)
CADop_diffr = CADop_diffr(~removeMask, :);

end