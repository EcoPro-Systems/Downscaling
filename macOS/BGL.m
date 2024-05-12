function BGL(path)
  % Import data from .csv files
  
  % Load observational data from Obs_data.csv
  Obs = table2array(readtable(strcat(path,'Obs_data.csv')));
  Obs_pixels = Obs(:, 1:2);
  Obs(:, 1:2) = [];

  % Load stand downscaled training data from Stand_Downscaled_training.csv
  Stand_SD = table2array(readtable(strcat(path, 'Stand_Downscaled_training.csv')));
  column_names=Stand_SD(1,:);
  Stand_SD(:, 1:2) = [];
  Stand_SD(1, :) = [];
  % Calculate dimensions
  N = size(Obs, 1);
  T_o = size(Obs, 2);

  % Load Mu1 and Mu2 data
  Mu1 = table2array(readtable(strcat(path, 'Mu1.csv')));
  Mu2 = table2array(readtable(strcat(path, 'Mu2.csv')));

  % Subtract the mean from the processes
  Z1 = Stand_SD(:, 1:T_o) - Mu1(:, 1:T_o);
  Z2 = Obs - Mu2(:, 1:T_o);

  % Remove unnecessary columns
  Stand_SD(:, 1:T_o) = [];
  Mu1(:, 1:T_o) = [];
  Mu2(:, 1:T_o) = [];
  
  % Calculate Z1 downscale and Z2 mean downscale
  Z1_downscale = Stand_SD - Mu1;
  Z2_mean_downscale = Mu2;

  % Define some constants and dimensions
  t_downscale = width(Z1_downscale);
  N = size(Obs, 1);

  % Setup Data

  % Create a data_array with two variables, N locations, and T_o realizations
  data_array = NaN([2, N, T_o]);
  data_array(1, :, :) = Z1;
  data_array(2, :, :) = Z2;
  p = size(data_array, 1);
  n = size(data_array, 2);
  m = size(data_array, 3);

  %% Basis Setup

  % Calculate orthogonal basis matrix using EOFs
  reshaped_data_array = horzcat(Z1, Z2);
  % Choose a number of basis functions to explain 99.99% of variance
  [U, S, ~] = svd(reshaped_data_array, 'econ');
  Cum_variance = NaN([width(S), 1]);
  for j = 1:width(S)
    S_diags = diag(S).^2 / sum(diag(S).^2);
    Cum_variance(j, 1) = sum(S_diags(1:j));
  end
  Var99_explained = find(Cum_variance >= 0.9999);
  U = U(:, 1:Var99_explained);

  % Retain basis functions for the stochastic term
  smallPhi = U;
  l = size(smallPhi, 2);

  %% Nugget Estimate

  % Calculate nugget estimates and smoothness estimates for each variable
  x0 = [1, 1];
  lb = [0, 0];
  ub = [50, 50];
  nuggetestimates = zeros([p, 1]);
  smoothnessestimates = zeros([p, 1]);

  for i = 1:p
    tt = smallPhi' * squeeze(data_array(i, :, :));
    smallPhi_S_Phi = tt * tt' / m;
    trS = norm(squeeze(data_array(i, :, :)), 'fro')^2 / m;
    [xval, fval] = fmincon(@(x) nuggetlikelihood_orthog(x, smallPhi_S_Phi, trS, n), x0, [], [], [], [], lb, ub);
    nuggetestimates(i) = xval(1);
    smoothnessestimates(i) = xval(2);
  end

  inv_error_variances = 1./nuggetestimates;
  
  %% Projecting Data and Creating Matrix

  % Project the data using the calculated nugget estimates
  projected_data = zeros(p * l, m);

  for i = 1:m
    projected_data(:, i) = vec((inv_error_variances .* data_array(:, :, i)) * smallPhi);
  end

  %% Main Algorithms

  % Calculate the BGL block diagonal matrix
  mleguess = BGLblockdiag_nopenalty(inv_error_variances, projected_data);

  % Calculate linear transformation matrices A1 and Q
  A1_start = eye(l);
  A1 = NaN([l, p * l]);
  for r = 1:l
    col1 = 2 * r - 1;
    col2 = 2 * r;
    A1(1:l, col1:col2) = horzcat(zeros(l, 1), A1_start(1:l, r));
  end
  Q = num2cell(mleguess, [1, 2]);
  Q = blkdiag(Q{:});

  % Calculate Q_inv
  mleguess_inv = NaN([p, p, l]);
  for r = 1:l
    mleguess_inv(:, :, r) = inv(mleguess(:, :, r));
  end
  Q_inv = num2cell(mleguess_inv, [1, 2]);
  Q_inv = blkdiag(Q_inv{:});

  Predicted_temp_all = NaN([N, t_downscale + 2]);
  Predicted_sd_all = NaN([N, t_downscale + 2]);

  Predicted_temp_all(:, 1:2) = Obs_pixels;
  Predicted_sd_all(:, 1:2) = Obs_pixels;

  % Calculate C and Sigma for future prediction
  for r = 1:l
    C(r, r) = mleguess_inv(1, 2, r) / mleguess_inv(1, 1, r);
    Sigma(r, r) = mleguess_inv(2, 2, r) - (mleguess_inv(1, 2, r))^2 / mleguess_inv(1, 1, r);
  end

  % Begin downscaling
  Start = cputime;
  tic
  for month = 1:t_downscale

    % New Z1 and the corresponding mean for future prediction
    Z1_new = Z1_downscale(:, month);
    Mean_new = Z2_mean_downscale(:, month);

    SC = smallPhi * C;
    Inv = inv(inv(A1 * Q_inv * A1') + eye(l) / nuggetestimates(1));

    SY = smallPhi' * Z1_new;
    Predicted_temp = Mean_new(:, 1) + (SC * Inv * SY) / nuggetestimates(1);

    Matrix = (Sigma + C * Inv * C');
    for loc = 1:n
      Predicted_sd = sqrt(smallPhi(loc, :) * Matrix * smallPhi(loc, :)' + nuggetestimates(2));
    end

    Predicted_temp_all(:, month + 2) = Predicted_temp;
    Predicted_sd_all(:, month + 2) = Predicted_sd;

    fprintf(strcat('t', num2str(month), ' done!'));
    fprintf('\n')
  end
  toc
  fprintf('\n')

  column_names(:,[1:2,(1:T_o)+2])=[];
  column_names=["lon","lat",column_names];
  Predicted_temp_all=array2table(Predicted_temp_all, 'VariableNames', column_names);
  Predicted_sd_all=array2table(Predicted_sd_all, 'VariableNames', column_names);
  
  writetable(Predicted_temp_all, strcat(path, '/BGL_Downscaled.csv'));
  writetable(Predicted_sd_all, strcat(path, '/BGL_Downscaled_UQ.csv'));

  % Delete unnecessary CSV files
  delete(strcat(path, '/Stand_Downscaled_training.csv'));
  delete(strcat(path, '/Mu1.csv'));
  delete(strcat(path, '/Mu2.csv'));
end
