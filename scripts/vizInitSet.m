function vizInitSet(bounds, eta)
    vertices = zeros(8, 3);
    vertices(1, :) = [bounds(1), bounds(3), eta];
    vertices(2, :) = [bounds(2), bounds(3), eta];
    vertices(3, :) = [bounds(2), bounds(4), eta];
    vertices(4, :) = [bounds(1), bounds(4), eta];
    vertices(5, :) = [bounds(1), bounds(3), 1];
    vertices(6, :) = [bounds(2), bounds(3), 1];
    vertices(7, :) = [bounds(2), bounds(4), 1];
    vertices(8, :) = [bounds(1), bounds(4), 1];

    faces = zeros(6, 4);
    faces(1, :) = [1, 2, 3, 4];
    faces(2, :) = [5, 6, 7, 8];
    faces(3, :) = [1, 2, 6, 5];
    faces(4, :) = [2, 3, 7, 6];
    faces(5, :) = [3, 4, 8, 7];
    faces(6, :) = [4, 1, 5, 8];

    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'cyan')
end