close all
clear
clc

lambdaArray = [
	0.9,  0.6; ...
	1.2,  0.8; ...
	1.2,  1.3; ...
	0.6,  -0.5; ...
	-0.6, -0.8; ...
	0.7 + 0.2 * 1i, 0.7 - 0.2 * 1i; ...
	-0.7 + 0.2 * 1i, -0.7 - 0.2 * 1i; ...
	];

for lambdaIterator = 1:size(lambdaArray, 1)
	lambdaIterator
	%% Tworzenie macierzy o zadanym widmie
	lambda = lambdaArray(lambdaIterator, :)
	A = createMatrix(lambda)

	%% Konstruowanie zbioru X0 punktów pocz¹tkowych po³o¿onych na okrêgu
	numberOfPoints = 16;
	r = 1;
	fi = 0:(2*pi)/numberOfPoints:2*pi - (2*pi)/numberOfPoints;
	X0 = [r * cos(fi); r * sin(fi)];

	%% Obliczanie trajektorii dla punktów pocz¹tkowych
	numberOfIterations = 7;
	trajectories = cell(numberOfPoints, 1);
	for i = 1:numberOfPoints
		trajectories{i} = calculateTrajectory(A, X0(:, i), numberOfIterations);
	end

	%% Portret fazowy
	figure(1)
		hold on
		scatter(X0(1, :), X0(2, :), 'b')
		for i=1:numberOfPoints
			plot(trajectories{i}(1, :), trajectories{i}(2, :), 'b')
			scatter(trajectories{i}(1, end), trajectories{i}(2, end), 'b*')
		end
		title("Trajektorie w przestrzeni stanów")
		legend('Punkty startowe', 'Trajektorie', 'Punkty koñcowe')
		axis equal
		pbaspect([1 1 1])
		trajectoryLimit = get(gca, 'xlim')
		
	%% Obraz macierzy A (A*X0)
	Y = A * X0;
	figure(2)
		hold on
		scatter(X0(1, :), X0(2, :), 'b', 'filled')
		scatter(Y(1, :), Y(2, :), 'r', 'filled') 
		title("Obraz Y = AX_0")
		legend("X_0", "Y")
		pbaspect([1 1 1])

	%% Ilustracja wartoœci i wektorów w³asnych
	[V, L] = eig(A);
	eigenValues = diag(L)
	eigenVectors = V

	if isreal(eigenValues)
		figure(3)
			hold on
			plot([0 eigenVectors(1, 1)], [0 eigenVectors(2, 1)], 'b:', 'LineWidth', 2)
			plot([0 eigenVectors(1, 2)], [0 eigenVectors(2, 2)], 'b:', 'LineWidth', 2)
			plot([0 eigenVectors(1, 1) * eigenValues(1)], [0 eigenVectors(2, 1) * eigenValues(1)], 'r', 'LineWidth', 1.25)
			plot([0 eigenVectors(1, 2) * eigenValues(2)], [0 eigenVectors(2, 2) * eigenValues(2)], 'r', 'LineWidth', 1.25)

			title("Ilustracja wartoœci i wektorów w³asnych")
			legend('v_1', 'v_2', '\lambda_1v_1', '\lambda_2v_2');
			pbaspect([1 1 1])
	end
	
	%% Ilustracja pola wektorowego
	samples = 2 * 10 + 1;
	linearSpace = linspace(floor(trajectoryLimit(1)), ceil(trajectoryLimit(2)), samples);
	[x, y] = meshgrid(linearSpace, linearSpace);
	startPoints = [x(:), y(:)]';
	endPoints   = (A * startPoints) - startPoints;
	xEnd = reshape(endPoints(1, :), samples, samples);
	yEnd = reshape(endPoints(2, :), samples, samples);

	figure(4)
		quiver(x, y, xEnd, yEnd, 'b')
		title("Pole wektorowe")
		xlim(trajectoryLimit)
		ylim(trajectoryLimit)
		pbaspect([1 1 1])

	%% Punkty pocz¹tkowe i koñcowe i wektory w³asne
	if isreal(eigenValues)
		numberOfPoints1 = 100;
		r = 1;
		fi = 0:(2*pi)/numberOfPoints1:2*pi - (2*pi)/numberOfPoints1;
		X1 = [r * cos(fi); r * sin(fi)];
		Y1 = A * X1;
		figure(5)
			hold on
			plot(X1(1, :), X1(2, :), 'b.')
			plot(Y1(1, :), Y1(2, :), 'r.') 

			plot([0 eigenVectors(1, 1)], [0 eigenVectors(2, 1)], 'b:', 'LineWidth', 2)
			plot([0 eigenVectors(1, 2)], [0 eigenVectors(2, 2)], 'b:', 'LineWidth', 2)
			plot([0 eigenVectors(1, 1) * eigenValues(1)], [0 eigenVectors(2, 1) * eigenValues(1)], 'r', 'LineWidth', 1.25)
			plot([0 eigenVectors(1, 2) * eigenValues(2)], [0 eigenVectors(2, 2) * eigenValues(2)], 'r', 'LineWidth', 1.25)

			title("Obraz Y = AX_0 oraz wektory i wartoœci w³asne")
			legend("X_0", "Y", 'v_1', 'v_2', '\lambda_1v_1', '\lambda_2v_2', 'FontSize', 9)
			
			pbaspect([1 1 1])
			axis equal
			imageLimit = get(gca, 'xlim')
			
			figure(2)
			xlim(imageLimit)
			ylim(imageLimit)
			figure(3)
			xlim(imageLimit)
			ylim(imageLimit)
	end
	
	%% Trajektorie na tle pola wektorowego
	figure(6)
		hold on
		scatter(X0(1, :), X0(2, :), 'b')
		for i=1:numberOfPoints
			plot(trajectories{i}(1, :), trajectories{i}(2, :), 'b')
			scatter(trajectories{i}(1, end), trajectories{i}(2, end), 'b*')
		end
		quiver(x, y, xEnd, yEnd, 'r')
		title("Trajektorie na tle pola wektorowego")
		legend('Punkty startowe', 'Trajektorie', 'Punkty koñcowe')
		pbaspect([1 1 1])
		xlim(trajectoryLimit)
		ylim(trajectoryLimit)
	
	%% Zapis figur
	for i = 1:6
		figure(i)
		saveas(gcf, ['fig/' num2str(lambdaIterator) '_' num2str(i) '.emf'], 'meta')
	end
	close all
end