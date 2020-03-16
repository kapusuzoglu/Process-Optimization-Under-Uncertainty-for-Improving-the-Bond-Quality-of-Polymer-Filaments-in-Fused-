import pymc3 as pm
import matplotlib.pyplot as plt
import seaborn as sns
import AM_Paper2_MCMC as fun1
# import numpy as np
# import covadapt

# id of the final part that wants to be included in the analysis
num_parts_last = 15

# Percent of the whole data is used to train the model
ratio = .05

# Call the measurements_and_training_data function to get
# x: Training data [printer nozzle temperature (t_n), speed (v_p), layer heigth (hth), x, y (coordinates of bonds)]
# y: bond length measurements at each interface of multiple parts
# numTotalParts: number of total parts analyzed (+1 is the real value because of indexing starts at 0)
# time_frames: time frame data of temperature measurements
# num_interfaces_of_each_part: numpy array stores number of interfaces of each part
# int_temperatures: interface temperature data of each part
# train_idx: randomly selected proportion (ratio) of training data
(X, y, numTotalParts, time_frames, int_temperatures, train_idx) = \
    fun1.measurements_and_training_data(num_parts_last, ratio)

print("Total number of interfaces (training data size): ", train_idx.shape[0], train_idx)

# Material parameters for the bond length model
surf_tension = 0.029  # Surface tension
b1_ = 3.45e-3  # Model parameter for temp. dependent surface tension
viscosity = 5100  # viscosity at 240 C

# Input parameters to the neck growth model
ModelParams = [surf_tension, b1_, viscosity]

# Parameters to be Calibrated
# 1- beta: neck growth model parameter
# 2- a0, a1, a2, a3: trend function coefficients (for printer temperature, speed, layer height, T_N, v_p, h)
# 3- length scale parameters of GP for 5 inputs, T_N, v_p, h, x & y coordinates
# 4- SigmaF - variance of GP
# In total, 11 parameters to calibrate.

# niter = 5000 # num. of chains

if __name__ == '__main__':
    with pm.Model() as model:

        # linear trend function coefficients:
        #  a0 + a1*x1 + a2*x2 + a3*x3 + a4*x4 + a5*x5, x = [T_N, v_p, hth, x, y]
        a0 = pm.Normal('a0', mu=0, sd=1)
        a = pm.Normal('a1', mu=0, sd=1, shape=5)

        # Priors for unknown neck growth model parameter
        # Normal distribution while constraining the support to be positive
        beta = pm.Bound(pm.Normal, lower=0.0)('beta', mu=0.056, sigma=0.05)

        # Observations: The low fidelity 1-D analytical model
        # Predicted_BL: bond length between two filaments of all parts:
        Predicted_BL = fun1.bond_model(beta, ModelParams, numTotalParts,
                                  int_temperatures, time_frames, train_idx)

        lin_trend_func = pm.gp.mean.Linear(coeffs=a, intercept=(Predicted_BL + a0))

        # GP length scale parameters for x = [T_N, v_p, hth, x, y] & process noise
        # fairly uninformative priors for the scale hyperparameters of the covariance functions,
        # and informative Gamma parameters for lengthscales
        l = pm.Gamma('l', alpha=1, beta=1, shape=5)
        # l = pm.Uniform('l', lower=1e-5, upper=10, shape=5)

        # sigma = pm.Uniform('sigma', lower=0.01, upper=3)
        eta = pm.HalfNormal('eta', sigma=2)

        # # uninformative prior on the drift amplitude
        # log_s2_d = pm.Uniform('log_s2_d', lower=-10, upper=5)
        # s2_d = pm.Deterministic('s2_d', T.exp(log_s2_d))
        #
        # # uninformative prior on the white noise variance
        # log_s2_w = pm.Uniform('log_s2_w', lower=-10, upper=5)
        # s2_w = pm.Deterministic('s2_w', T.exp(log_s2_w))

        # Specify the covariance functions for each Xi
        # Since the covariance is a product, only scale one of them by eta.
        # Scaling both overparameterizes the covariance function.
        cov = eta ** 2 * pm.gp.cov.ExpQuad(input_dim=5, ls=l)

        # The scale of the white noise term can be provided,
        sigma = pm.HalfNormal('sigma', sigma=2)

        # a_non = pm.Deterministic('a_non', )


        gp = pm.gp.Marginal(mean_func=lin_trend_func, cov_func=cov)
        y_ = gp.marginal_likelihood('y_', X=X, y=y, noise=sigma)

        # inference
        # vi_fit = pm.fit(method='advi', n=100000)
        # trace = vi_fit.sample(100000)

        # trace = pm.sample(5000, init='advi', n_init=50000)

        # # Define the potential
        # pot = covadapt.EigvalsAdapt(
        #     model.ndim,
        #     np.zeros(model.ndim),
        #     estimators=[
        #         lambda samples, grads:
        #         covadapt.eigh_lw_samples_grads(
        #             samples, grads, n_eigs=20, n_eigs_grad=20, n_final=40
        #         )
        #     ],
        #     display=True,
        # )
        #
        # # Initialize the NUTS sampler with the new potential
        # step = pm.NUTS(potential=pot)

        step = pm.NUTS()  # Hamiltonian MCMC with No U-Turn Sampler
        # step = pm.Metropolis()
        trace = pm.sample(40000, step, tune=1000, cores=1, chains=2)
                          # progressbar=True)
        # trace = pm.sample(20000, step, tune=1000, cores=1, chains=2,
        #                   init='advi+adapt_diag',
        #                   progressbar=True)

        pm.traceplot(trace)
        plt.show()

        prior = pm.sample_prior_predictive(10000)

        var_names = [j for j in pm.util.get_default_varnames(trace.varnames, include_transformed=False) if trace[j].ndim
                     == 1]


fig, axes = plt.subplots(len(var_names), 1, figsize=(12, 16))
for ax, var_name in zip(axes.flatten(), var_names):
    for label, obj in (('posterior', trace), ('prior', prior)):
        sns.distplot(obj[var_name], hist=False, kde=True, ax=ax, label=label, norm_hist=True)
    ax.legend()
    ax.set_title(var_name)
plt.show()


        # ax = plt.tight_layout()
        # format_axes(ax)
        # plt.savefig("../image.pdf")

    # pm.summary(trace).round(2)

    # Calculate effective sample size: \hat{n}_{eff} = \frac{mn}{1 + 2 \sum_{t=1}^T \hat{\rho}_t}
    # pm.effective_n(trace)

    # Evaluate convergence
    # Gelman-Rubin: \hat{R} = \frac{\hat{V}}{W}

# Open file and write
# import _pickle as pickle
# file_temp = open('traces20000', 'wb')
# pickle.dump(trace, file_temp)
# file_temp.close()
#
# Open file and read
# file_temp = open('traces20000', 'rb')
# traces = pickle.load(file_temp)
# file_temp.close()


# end = time.time()
# print(end-start, "seconds")

# detailed summary of the posterior distribution for each parameter
# pm.plot_posterior(trace);

# # We can also consider all the variables simultaneously.
# pm.autocorrplot(trace, max_lag=20);

# # From these plots we see that the auto-correlation is not problematic. Indeed, we can test this through the
# # effective sample size, which should be close to the total number of samples
# pm.diagnostics.effective_n(trace)

# plt.plot(x, posterior(x,y))
# plt.plot(x, ss.gamma.pdf(x,a=a,scale=1/b), 'r-')
# plt.title('Prior and Posterior Distributions')
# plt.show()
