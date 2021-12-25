function [bootstat] = matrixBootstrap(nboot, bootfun, mat)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

n_obsv = size(mat, 1);
n_points = size(mat, 2);

perm = randi(n_obsv, nboot, n_points);

bootmat = zeros(nboot, n_points);

for i=1:nboot
    for j=1:n_points
        bootmat(i, j) = mat(perm(i, j), j);
    end
end

bootstat = feval(bootfun, bootmat, 2);

end

