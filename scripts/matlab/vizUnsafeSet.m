function vizUnsafeSet(bounds)
    cap = 5.5;
    vertices = zeros(8, 3);
    vertices(1, :) = [bounds(1), bounds(3), 0];
    vertices(2, :) = [bounds(2), bounds(3), 0];
    vertices(3, :) = [bounds(2), bounds(4), 0];
    vertices(4, :) = [bounds(1), bounds(4), 0];
    vertices(5, :) = [bounds(1), bounds(3), .1];
    vertices(6, :) = [bounds(2), bounds(3), .1];
    vertices(7, :) = [bounds(2), bounds(4), .1];
    vertices(8, :) = [bounds(1), bounds(4), .1];
    vertices(9, :) = [bounds(1), bounds(3), cap];
    vertices(10, :) = [bounds(2), bounds(3), cap];
    vertices(11, :) = [bounds(2), bounds(4), cap];
    vertices(12, :) = [bounds(1), bounds(4), cap];

    faces_sld = zeros(6, 4);
    faces_sld(1, :) = [1, 2, 3, 4];
    faces_sld(2, :) = [5, 6, 7, 8];
    faces_sld(3, :) = [1, 2, 6, 5];
    faces_sld(4, :) = [2, 3, 7, 6];
    faces_sld(5, :) = [3, 4, 8, 7];
    faces_sld(6, :) = [4, 1, 5, 8];

    faces_op = zeros(6, 4);
    faces_op(1, :) = [1, 2, 3, 4] + 4;
    faces_op(2, :) = [5, 6, 7, 8] + 4;
    faces_op(3, :) = [1, 2, 6, 5] + 4;
    faces_op(4, :) = [2, 3, 7, 6] + 4;
    faces_op(5, :) = [3, 4, 8, 7] + 4;
    faces_op(6, :) = [4, 1, 5, 8] + 4;

    patch('Vertices', vertices, 'Faces', faces_sld, 'FaceColor', 'red')
    % patch('Vertices', vertices, 'Faces', faces_op, 'FaceColor', 'red', 'FaceAlpha', 0.2)
end

