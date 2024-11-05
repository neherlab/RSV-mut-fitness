import math
import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize
from matplotlib import pyplot as plt


def get_probabilistic_estimate(
    n_obs,
    n_pred,
    tau_squared,
    show_p_f_given_n_obs_n_pred=False,
    offset=1.0e-6,  # logarithm offset
    s_prior=2.5,  # gaussian prior std dev
    N_f=150,  # number of points for fitness grid
    N_n_exp=300,  # number of points for expected counts grid
):
    # Gaussian prior variance
    s2_prior = s_prior * s_prior

    # P(f) (normalized)
    def p_f(f):
        # Gaussian with mean 0 and standard deviation 2.5
        return norm.pdf(f, loc=0.0, scale=s_prior)

    # Precompute constants used in P(n_exp|n_pred)
    c2 = 1 / math.sqrt(tau_squared * 2 * math.pi)
    c3 = 2 * tau_squared

    # P(n_exp|n_pred) (normalized)
    def p_n_exp_given_n_pred(n_exp):
        # Log-normal distribution with mean log(n_pred) and variance tau^2
        return (
            c2
            * np.exp(-(np.log((n_exp + offset) / n_pred) ** 2) / c3)
            / (n_exp + offset)
        )

    # Precompute constant used in P(n_obs|f,n_exp)
    log_n_obs_factorial = math.log(math.factorial(n_obs))

    # P(n_obs|f,n_exp) (normalized)
    def p_n_obs_given_f_n_exp(n_exp, f):
        # Poisson distribution with mean n_exp * e^f
        log_p = (
            n_obs * (np.log(n_exp + offset) + f)
            - n_exp * np.exp(f)
            - log_n_obs_factorial
        )
        return np.exp(log_p)

    # Precompute constants used in P(n_obs|f,n_pred)
    mean_log_normal = np.exp(np.log(n_pred) + tau_squared / 2)
    std_log_normal = np.sqrt(
        (np.exp(tau_squared) - 1) * np.exp(2 * np.log(n_pred) + tau_squared)
    )
    upper_bound_log_normal = mean_log_normal + 10 * std_log_normal
    lower_bound_log_normal = mean_log_normal - 10 * std_log_normal
    n_obs_plus_1 = 1 + n_obs
    n_obs_plus_1_sqrt = np.sqrt(1 + n_obs)

    # P(n_obs|f,n_pred) = int_{0}^{inf} d{n_exp} P(n_obs|f,n_exp)P(n_exp) (normalized)
    def p_n_obs_given_f_n_pred(f):
        # Set the numerical integration boundaries according to the shape of the involved distributions
        e_to_the_minus_f = np.exp(-f)
        mean_poisson = e_to_the_minus_f * n_obs_plus_1
        std_poisson = e_to_the_minus_f * n_obs_plus_1_sqrt
        lower_lim = max(0, max(lower_bound_log_normal, mean_poisson - 5 * std_poisson))
        upper_lim = min(upper_bound_log_normal, mean_poisson + 5 * std_poisson)

        # Plot integrand
        # plot_integrand(f, lower_lim, upper_lim)

        # integrand: P(n_exp|n_pred)*P(n_obs|f,n_exp)
        n_exp_values = np.linspace(lower_lim, upper_lim, 100).flatten()
        integrand_values = p_n_exp_given_n_pred(n_exp_values) * p_n_obs_given_f_n_exp(
            n_exp_values, f
        )

        # Perform the integration over n_exp using the trapezoidal rule
        integral = np.trapz(y=integrand_values, x=n_exp_values)

        return integral

    # P(f|n_obs,n_pred) ~ P(n_obs|f,n_pred)*P(f) (normalization 1/P(n_obs|n_pred) omitted)
    def p_f_given_n_obs_n_pred(f):
        return p_n_obs_given_f_n_pred(f) * p_f(f)

    # f_hat = argmax_f P(f|n_obs,n_pred)
    def get_f_hat():
        # Determine f_hat by finding the argmin of -log[P(f|n_obs,n_pred)]. Use the naive estimate as an initial guess.
        initial_guess = np.log(n_obs + 0.5) - np.log(n_pred + 0.5)
        result = minimize(
            lambda f: -np.log(p_f_given_n_obs_n_pred(f) + offset),
            x0=initial_guess,
            method="BFGS",
        )

        return result.x[0]

    # Plot P(n_exp)*P(n_obs|f,n_exp) as a function of n_exp
    def plot_integrand(f, lower_lim, upper_lim):
        n_exp_values = np.linspace(lower_lim, upper_lim, 500)
        integrand_values = [
            p_n_exp_given_n_pred(n_exp) * p_n_obs_given_f_n_exp(n_exp, f)
            for n_exp in n_exp_values
        ]

        plt.figure(figsize=(8, 4.5))
        plt.plot(
            n_exp_values,
            integrand_values,
            label=r"Integrand $p(n_{exp}) \cdot p(n_{obs}|n_{exp}, f)$",
        )
        plt.xlabel("$n_{exp}$")
        plt.ylabel("Integrand")
        plt.title(f"f = {f}")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    # Compute mean, highest density interval, and standard deviation of P(f|n_obs,n_pred)
    def characterize_p_f_given_n_obs_n_pred(f_max):
        # Define resolution of integration domain
        # N_f, N_n_exp = 150, 300

        # Define f integration domain (using the analytical approximation of the posterior)
        mu_posterior = (
            np.log((n_obs + offset) / n_pred)
            * s2_prior
            / (tau_squared + 1 / (n_obs + offset) + s2_prior)
        )
        sigma_posterior = np.sqrt(
            1 / (1 / (tau_squared + 1 / (n_obs + offset)) + 1 / s2_prior)
        )
        f_lower_lim, f_upper_lim = (
            mu_posterior - 5 * sigma_posterior,
            mu_posterior + 5 * sigma_posterior,
        )
        f_values = np.linspace(f_lower_lim, f_upper_lim, N_f)

        # Define n_exp integration domain (using the shape of the involved distributions)
        e_to_the_minus_f_lower = np.exp(-f_lower_lim)
        e_to_the_minus_f_upper = np.exp(-f_upper_lim)
        mean_poisson_lower = e_to_the_minus_f_lower * n_obs_plus_1
        std_poisson_lower = e_to_the_minus_f_lower * n_obs_plus_1_sqrt
        mean_poisson_upper = e_to_the_minus_f_upper * n_obs_plus_1
        std_poisson_upper = e_to_the_minus_f_upper * n_obs_plus_1_sqrt
        n_exp_lower_lim = max(
            0,
            lower_bound_log_normal,
            min(
                mean_poisson_lower - 5 * std_poisson_lower,
                mean_poisson_upper - 5 * std_poisson_upper,
            ),
        )
        n_exp_upper_lim = min(
            upper_bound_log_normal,
            max(
                mean_poisson_lower + 5 * std_poisson_lower,
                mean_poisson_upper + 5 * std_poisson_upper,
            ),
        )
        n_exp_values = np.linspace(n_exp_lower_lim, n_exp_upper_lim, N_n_exp)

        # Create 2D mesh (x: f, y: n_exp) over integration domain
        n_exp_mesh, f_mesh = np.meshgrid(n_exp_values, f_values, indexing="ij")

        # Compute P(n_exp|n_pred)*P(n_obs|f,n_exp)*P(f) on the entire mesh
        integrand_values = (
            p_n_exp_given_n_pred(n_exp_mesh)
            * p_n_obs_given_f_n_exp(n_exp_mesh, f_mesh)
            * p_f(f_mesh)
        )

        # Integrate over n_exp
        n_exp_integral = np.trapz(integrand_values, x=n_exp_values, axis=0)

        # Integrate over f to get the normalization constant
        p_n_obs_given_n_pred = np.trapz(n_exp_integral, x=f_values)

        # Calculate the mean of the posterior
        integrand_values = f_values * n_exp_integral
        f_mean = np.trapz(integrand_values, x=f_values) / p_n_obs_given_n_pred

        # Calculate the standard deviation of the posterior
        integrand_values = (f_values - f_mean) ** 2 * n_exp_integral
        f_std = np.sqrt(np.trapz(integrand_values, x=f_values) / p_n_obs_given_n_pred)

        # Get the highest density interval of the posterior
        left_interval, right_interval = get_hdi(
            f_values, n_exp_integral / p_n_obs_given_n_pred
        )

        # Plot p_f_given_n_obs(f) to validate result
        if show_p_f_given_n_obs_n_pred:
            plt.figure(figsize=(8, 4.5))
            plt.plot(
                f_values, n_exp_integral / p_n_obs_given_n_pred, label=r"$p(f|n_{obs})$"
            )
            plt.axvline(x=f_hat, color="green", label=r"new $\hat{f}$")
            # plt.axvline(x=f_mean, color='blue', label=r'mean')
            orig_estim = np.log((n_obs + 0.5) / (n_pred + 0.5))
            plt.axvline(x=orig_estim, color="red", label=r"old $\hat{f}$")
            plt.fill_between(
                f_values,
                n_exp_integral / p_n_obs_given_n_pred,
                where=((f_values >= left_interval) & (f_values <= right_interval)),
                color="gray",
                alpha=0.5,
            )
            plt.xlabel("$f$")
            plt.legend()
            plt.grid(True)
            plt.title(
                rf"n_pred = {round(n_pred, 2)}, n_obs = {round(n_obs, 2)}, $\tau^{2}$ = {round(tau_squared, 2)}"
            )
            plt.tight_layout()
            plt.show()

        return f_mean, f_max - left_interval, right_interval - f_max, f_std

    # Find f_hat = argmax_f P(f|n_obs,n_pred)
    f_hat = get_f_hat()

    # Get mean, highest density interval, and standard deviation of P(f|n_obs,n_pred)
    f_mean, left_confidence_interval, right_confidence_interval, f_std = (
        characterize_p_f_given_n_obs_n_pred(f_hat)
    )

    return f_hat, left_confidence_interval, right_confidence_interval, f_mean, f_std


def get_hdi(x, p, coverage=0.682):
    # Calculate the contribution of all intervals
    dx = x[1] - x[0]
    dp = dx * (p[:-1] + (p[1:] - p[:-1]) / 2)

    # Sort x and corresponding dp in descending order of dp
    sorted_indices = np.argsort(-dp)
    x_sorted = x[sorted_indices]
    dp_sorted = dp[sorted_indices]

    # Compute cumulative sum of probabilities
    cumulative_prob = np.cumsum(dp_sorted)

    # Find the smallest set of points that cover the desired probability
    cutoff_idx = np.argmax(cumulative_prob >= coverage)
    hdi_x_values = x_sorted[: cutoff_idx + 1]

    # Get the smallest interval by taking the min and max of the selected x values
    interval = (min(hdi_x_values), max(hdi_x_values))

    return interval


def add_probabilistic_estimates(dataframe, offset=1.0e-6, N_f=150, N_n_exp=300):
    def calculate_estimates(row):
        return get_probabilistic_estimate(
            row.actual_count,
            row.predicted_count,
            row.tau_squared,
            show_p_f_given_n_obs_n_pred=False,
            offset=offset,
            N_f=N_f,
            N_n_exp=N_n_exp,
        )

    estimates = dataframe.apply(calculate_estimates, axis=1)

    (
        dataframe["f_max"],
        dataframe["left_conf_int"],
        dataframe["right_conf_int"],
        dataframe["f_mean"],
        dataframe["f_st_dev"],
    ) = zip(*estimates)
