function mat_axes = matrix_plot(data, qoi_text, plot_params)
% MATRIX_PLOT: Function generates a matrix plot of the parameter 2D
% contours for each dimension based on Kernel Density Estimation
%
% See COPYRIGHT for license information 
%

% Determine the number of dimensions and datapoints
data = data';
N_dims = size(data,1);

% Tolerance for variable bound
tol = 0;

% Resolution for plotting
N_points = 40;

%% Determine Bounds for Variables

% Declare matrix to store bounds
bounds = zeros(2,N_dims);

% Find all bounds
for i=1:N_dims
    
    % Find bounds
    temp_min = min(data(i,:));
    temp_max = max(data(i,:));
    offSet = tol*(temp_max - temp_min)/2;
    iymin  = temp_min - offSet;
    iymax  = temp_max + offSet;

    % Save in bounds
    bounds(1,i) = iymin;
    bounds(2,i) = iymax;

end

%% Generate all subplots

% Mesh for plotting
[x0 y0] = meshgrid(linspace(0,1,N_points),linspace(0,1,N_points));
x1 = reshape(x0,N_points*N_points,1);
y1 = reshape(y0,N_points*N_points,1);

% Right Subplot axes
Nh = N_dims;
Nw = N_dims;
gap    = [.01 .01];
marg_h = [.1 .1];
marg_w = [.1 .1];
axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
py = 1-marg_h(2)-axh; 

figure;
hold on

mat_axes = cell(N_dims, N_dims);

for i=1:N_dims
    
    px = marg_w(1);
    for j=1:i

        subplot('Position',[px py axw axh]);
        px = px+axw+gap(2);
        
        hold on
        
        % If j==i (then plot the density)
        if j==i
            
            kde_eval = kde(data(i,:),'rot');
            xs_temp  = linspace(bounds(1,i),bounds(2,i),2*N_points);
            ys_temp  = evaluate(kde_eval,xs_temp);
            plot(xs_temp,ys_temp);
            xlim([bounds(1,i) bounds(2,i)])
            title(['$y_{' num2str(i) '}$'],'interpreter','latex')

        else
            
            % Evaluate bounds
            x = x1.*(bounds(2,i)-bounds(1,i)) + bounds(1,i);
            y = y1.*(bounds(2,j)-bounds(1,j)) + bounds(1,j);
            
            % Evaluate KDE
            temp = kde(data([i j],:),'rot');
            
            % Generate contour
            contour(reshape(y, N_points, N_points), ...
                    reshape(x, N_points, N_points), ...
                    reshape(evaluate(temp, [x y]'), N_points, N_points))
            xlim([bounds(1,j) bounds(2,j)])
            ylim([bounds(1,i) bounds(2,i)])
            colormap winter

        end
        % Plot details
        if i==j
            title(qoi_text{i},'interpreter','latex')
        end

        % Remove ylabels if plot is on the interior of matrix
        if(j~=1)
            set(gca,'yticklabel',[])
        else
            set(gca,'TickLabelInterpreter','latex','FontWeight','normal')
        end
                
        % Remove xlabels if plot is on the interior of matrix
        if(i~=N_dims)
            set(gca,'xticklabel',[])
        else
            set(gca,'TickLabelInterpreter','latex','FontWeight','normal')
        end
        
        % Multiply Title Font Size
        set(gca,'TitleFontSizeMultiplier', 1.2, 'FontSize', plot_params.fsz, 'LineWidth', plot_params.alw);
        hold off

        mat_axes{i,j} = get(gca,'XLim'); %get(gca,'XLim'), get(gca,'YLim');
        
    end
    py = py-axh-gap(1);
end

hold off

end