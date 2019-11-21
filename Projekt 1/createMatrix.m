function A = createMatrix(lambda)
	A = eye(2);
	
	while true
		V = rand(2, 2) - 0.5;
		if abs(det(V)) > 0.1
			break
		end
	end
	
	if prod(size(lambda) == [1 2]) || prod(size(lambda) == [2 1])
		if isreal(lambda(1)) && isreal(lambda(2))
			A = diag(lambda);
			A = V * A / V;
		elseif conj(lambda(1)) == lambda(2)
			A = diag(real(lambda)) + flip(diag(imag(lambda)));
			A = V * A / V;
		end
	end
end

