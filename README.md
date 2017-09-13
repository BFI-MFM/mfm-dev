# MFM Model Repository

The MFM Model Repository is an ongoing project that aims to assess and include models and frameworks in the macro-financial literature. It also evaluates numerical and computational methods that are essential in solving models.

To download and install the latest version, use the follow command line:

<code>git clone https://github.com/BFI-MFM/mfm-dev.git</code>

You should see directory <code>mfm_dev</code> being created.  Run the following at the command line:

```
cd mfm-dev
git submodule init
git submodule update
```

After that, start MATLAB and run

```matlab
set_path_mfm
```

This will set up the environment for you.

## Models

Currently, we have a few models available ([Brunnermeier & Sannikov (2014)](/brusan), [He & Krishnamurthy (2013)](/hk), [Brock & Mirman (1972)](/brock_mirman), and [Klimenko, Pfeil, Rochet, & Nicolo (2016)](/rochet)). As an example, we will show how you can run Brunnermei and Sannikov (2014).

### Brunnermeier and Sannikov (2014)

The model is based on *A Macroeconomic Model with a Financial Sector*. A more detailed instruction file is in the [model folder](/brusan). To compute the equilibrium, run
```matlab
solve_equilibrium
```
The model will generate two figures:

* Figure 1: shows important equilibrium quantities such as q, psi, and drift and volatility of eta.

* Figure 2: shows expert and household utility within the model.

You can customize the model parameters. Please read the detailed instructions in the model folder.

## Tools and Numerical Methods

As of now we have developed tools to compute shock elasticities based on the work by Borovicka and Hansen (2012) and the numerical method to compute Chernoff entropy in Hansen and Sargent (2012). We will use shock elasticities as an example.

### Shock Elasticities

Both shock exposure and price elasticities can be computed, along with many auxiliary measures such as the stationary density and term structures, through the toolbox deveoped by Yiran Fan and Paymon Khorrami. As a demo, you can run

```matlab
demo_BS2014_SEimfd1
```

It uses the model by Brunnermeier & Sannikov (2014) as an example and computes the shock elasticities and other auxiliary outputs (stationary density, etc.)

The toolbox is versatile and contains many useful functions. For more information, look at the [PDF guide](/shockElas/main.pdf).

## To Do List

The MFM team is very active in assessing models and developing code based on the macro-financial literature. In the near future, we aim to add
* More macro-financial models in different programming languages
* A model comparison framework to diagnose and assess model strengths and weaknesses
* Tools to measure uncertainty and risk in models
* High performance computing tools for solving models
