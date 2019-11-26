close all
clear
clc

lambdaArray = [
	1, 1;
	0.6, 0.6;
	1.5, 1.5;
	-0.5, -0.5;
	0, 0;
	
	0.7,  0.6;
	1.2,  0.8;
	1.2,  1.3;
	
	0.6,  -0.5;
	0.5,  -1.2;
	1.2,  -0.4;
	-0.6, -0.8;
	-0.6, -1.2;
	-1.4, -1.2;
	
	0.5,  1;
	1.5,  1;
	-0.5, 1;
	-1.5, 1;
	
	0.8,  0;
	1.3,  0;
	-0.5, 0;
	-1.2, 0;
	
	0.5 + 0.3 * 1i, 0.5 - 0.3 * 1i;
	1.2 + 0.2 * 1i, 1.2 - 0.2 * 1i;
	-0.5 + 0.4 * 1i, -0.5 - 0.4 * 1i;
	-1.3 + 0.4 * 1i, -1.3 - 0.4 * 1i;
	0 + 0.5 * 1i, 0 - 0.5 * 1i;
	0 + 1 * 1i, 0 - 1 * 1i;
	0.9 + 0.9 * 1i, 0.9 - 0.9 * 1i;
	];

for lambdaIterator = 23:size(lambdaArray, 1)
	lambdaIterator;
	%% Tworzenie macierzy o zadanym widmie
	lambda = lambdaArray(lambdaIterator, :)
	A = createMatrix(lambda);

	%% Konstruowanie zbioru X0 punktów pocz¹tkowych po³o¿onych na okrêgu
	numberOfPoints = 16;
	r = 1;
	fi = 0:(2*pi)/numberOfPoints:2*pi - (2*pi)/numberOfPoints;
	X0 = [r * cos(fi); r * sin(fi)];

	%% Obliczanie trajektorii dla punktów pocz¹tkowych
	numberOfIterations = 4;
	trajectories = cell(numberOfPoints, 1);
	for i = 1:numberOfPoints
		trajectories{i} = calculateTrajectory(A, X0(:, i), numberOfIterations);
	end

	%% Portret fazowy
	fig = figure;
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
		trajectoryLimit = get(gca, 'xlim');
		close(fig)
		
	%% Obraz macierzy A (A*X0)
	Y = A * X0;

	%% Ilustracja wartoœci i wektorów w³asnych
	[V, L] = eig(A);
	eigenValues = diag(L)
	eigenVectors = V
	
	%% Ilustracja pola wektorowego
	samples = 2 * 5 + 1;
	linearSpace = linspace(floor(trajectoryLimit(1)), ceil(trajectoryLimit(2)), samples);
	[x, y] = meshgrid(linearSpace, linearSpace);
	startPoints = [x(:), y(:)]';
	endPoints   = (A * startPoints) - startPoints;
	xEnd = reshape(endPoints(1, :), samples, samples);
	yEnd = reshape(endPoints(2, :), samples, samples);

	%% Punkty pocz¹tkowe i koñcowe i wektory w³asne
	numberOfPoints1 = 100;
	r = 1;
	fi = 0:(2*pi)/numberOfPoints1:2*pi - (2*pi)/numberOfPoints1;
	X1 = [r * cos(fi); r * sin(fi)];
	Y1 = A * X1;
	figure
		hold on
		plot([X1(1, :) X1(1, 1)], [X1(2, :) X1(2, 1)], 'b', 'LineWidth', 1.5)
		plot([Y1(1, :) Y1(1, 1)], [Y1(2, :) Y1(2, 1)], 'r', 'LineWidth', 1.5) 
		if isreal(eigenValues)
			plot([0 eigenVectors(1, 1)], [0 eigenVectors(2, 1)], 'm:', 'LineWidth', 3)
			plot([0 eigenVectors(1, 2)], [0 eigenVectors(2, 2)], 'm:', 'LineWidth', 3)
			plot([0 eigenVectors(1, 1) * eigenValues(1)], [0 eigenVectors(2, 1) * eigenValues(1)], 'g', 'LineWidth', 2)
			plot([0 eigenVectors(1, 2) * eigenValues(2)], [0 eigenVectors(2, 2) * eigenValues(2)], 'g', 'LineWidth', 2)
		end

		for i=1:numberOfPoints
			plot(trajectories{i}(1, :), trajectories{i}(2, :), '-o')
		end
		quiver(x, y, xEnd, yEnd, 'Color', 0.5*[1 1 1])

		title({["\lambda_1 = " + num2str(lambda(1)) + ", \lambda_2 = " +  num2str(lambda(2))]  ...
			["v_1 = [" + num2str(eigenVectors(1, 1)) + ", " + num2str(eigenVectors(2, 1)) + "]"]  ...
			["v_2 = [" + num2str(eigenVectors(1, 2)) + ", " + num2str(eigenVectors(2, 2)) + "]"]})
		
		if isreal(eigenValues)
			legend("X_0", "Y", 'v_1', 'v_2', '\lambda_1v_1', '\lambda_2v_2', 'FontSize', 9, 'Location', 'NorthEastOutside')
		else
			legend("X_0", "Y", 'FontSize', 9, 'Location', 'NorthEastOutside')
		end
		pbaspect([1 1 1])
		axis equal
		imageLimit = get(gca, 'xlim');
		xlim(trajectoryLimit);
		ylim(trajectoryLimit);

		saveas(gcf, ['fig/' num2str(lambdaIterator) '_combined.emf'], 'meta')

end