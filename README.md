# MFM Model Repository

The MFM Model Repository is an ongoing project that aims to assess and include models and framework in the macro-financial literature. It also evaluates numerical and computational methods that are essential in solving models.

To download and install the latest version, use the follow command line:

<code>git clone https://github.com/BFI-MFM/mfm-dev.git</code>

After that, start MATLAB and run

```matlab
set_path_mfm
```
This will set up the environment for you.

## Running models

Currently, we have models from Brunnermeier & Sannikov (2014) and He & Krishnamurthy (2013) available.

### Brunnermeier and Sannikov (2014)

The model is based on *A Macroeconomic Model with a Financial Sector*. A more detailed instruction file is in the [model folder](../tree/master/brusan). To compute the equilibrium, run
```matlab
solve_equilibrium
```
The model will generate two figures:

⋅⋅* Figure 1: shows important equilibrium quantities such as q, $\psi$, and drift and volatility of $\eta$.

⋅⋅* Figure 2: shows expert and household utility within the model.

You can customize the model parameters. Please read the detailed instructions in the model folder.
