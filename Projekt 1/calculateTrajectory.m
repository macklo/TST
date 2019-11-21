function x = calculateTrajectory(A, x0, numberOfIterations)
	x = zeros(2, numberOfIterations + 1);
	x(:, 1) = x0;
	for i = 2:numberOfIterations + 1
		x(:, i) = A * x(:, i - 1);
	end
end

