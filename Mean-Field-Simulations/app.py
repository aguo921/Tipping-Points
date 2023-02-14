import streamlit as st
import matplotlib.pyplot as plt
from matplotlib import cm
from simulation_functions import *
from simulation_plotting import *

@st.cache
def ensemble_paths_cache(*args, **kwargs):
    return ensemble_paths(*args, **kwargs)

@st.cache
def escaped_trajectories_cache(*args):
    return escaped_trajectories(*args)

@st.cache
def distribution_cache(*args, **kwargs):
    return distribution(*args, **kwargs)

@st.cache
def estimate_variance_cache(*args, **kwargs):
    return estimate_variance(*args, **kwargs)

@st.cache
def analytical_variance_cache(*args, **kwargs):
    return analytical_variance(*args, **kwargs)

@st.cache
def probability_density_cache(*args, **kwargs):
    return probability_density(*args, **kwargs)


def plot_grid(types, xlabel, ylabel, f):
    num_types = len(types)

    fig, axs = plt.subplots(ncols=num_types)

    for i, type in enumerate(types):
        ax = axs[i] if num_types > 1 else axs

        f(ax, type)
    
        ax.set_title(type)
        ax.set_xlabel(xlabel)
        if i == 0:
            ax.set_ylabel(ylabel)

    fig.set_figheight(4.8 - num_types*0.8)
    fig.tight_layout()

    return fig


st.title("Simulations of Slow-Fast Dynamical Systems")

st.latex(r"dx=f(x,y)dt+\sigma dW")
st.latex(r"dy=\epsilon dt")

st.header("Simulation Settings")

noise = st.number_input("Sigma", min_value=0., value=0.1)  # noise level
scale = st.number_input("Epsilon", min_value=0., value=0.02)  # time scale separation

y0 = st.number_input("y0", max_value=0., value=-1., step=0.1)  # initial y
yf = st.number_input("yf", min_value=y0, value=0., step=0.1)  # final y

T = (yf - y0) / scale  # run duration

h = st.number_input("h", min_value=0., value=0.1)  # step size

# type of bifurcations
types = ["Fold", "Transcritical", "Pitchfork"]

# normal form function
f = {
    "Fold": lambda x, y: -y-x**2,
    "Transcritical": lambda x, y: y*x-x**2,
    "Pitchfork": lambda x, y: y*x+x**3
}

# initial x condition at equilibrium
x0 = {
    "Fold": np.sqrt(-y0),
    "Transcritical": 0,
    "Pitchfork": 0
}

# stable and unstable equilibria
equilibria = {
    "Fold": [
        {   # unstable equilibria before bifurcation point
            "range": (y0, 0),
            "func": lambda y: -np.sqrt(-y),
            "stable": False
        },
        {   # stable equilibria before bifurcation point
            "range": (y0, 0),
            "func": lambda y: np.sqrt(-y),
            "stable": True
        },
    ],
    "Transcritical": [
        {   # unstable equilibria before bifurcation point
            "range": (y0, 0),
            "func": lambda y: y,
            "stable": False
        },
        {   # stable equilibria before bifurcation point
            "range": (y0, 0),
            "func": lambda y: 0,
            "stable": True
        },
        {   # unstable equilibria after bifurcation point
            "range": (0, yf),
            "func": lambda y: y,
            "stable": True
        },
        {   # stable equilibria after bifurcation point
            "range": (0, yf),
            "func": lambda y: 0,
            "stable": False
        }
    ],
    "Pitchfork": [
        {   # unstable equilibria before bifurcation point
            "range": (y0, 0),
            "func": lambda y: -np.sqrt(-y),
            "stable": False
        },
        {   # unstable equilibria before bifurcation point
            "range": (y0, 0),
            "func": lambda y: np.sqrt(-y),
            "stable": False
        },
        {   # stable equilibria before bifurcation point
            "range": (y0, 0),
            "func": lambda y: 0,
            "stable": True
        }
    ]
}

# probability density function p(x, y)
pdf = {
    "Fold": lambda x, y: np.exp(2/noise**2*(-y*x - 1/3*x**3 + 2/3*(-y)**(3/2))),
    "Transcritical": lambda x, y: np.exp(2/noise**2*(1/2*y*x**2 - 1/3*x**3 - 1/6*y**3)),
    "Pitchfork": lambda x, y: np.exp(2/noise**2*(1/2*y*x**2 + 1/4*x**4 + 1/4*y**2))
}

# valid domain of x values for pdf given y value
domain = {
    "Fold": lambda y: (-np.sqrt(-y), 5),
    "Transcritical": lambda y: (y, 5),
    "Pitchfork": lambda y: (-np.sqrt(-y), np.sqrt(-y))
}

# x limits on graph
xlim = {
    "Fold": (-1, 1+3*noise),
    "Transcritical": (-1, 1),
    "Pitchfork": (-1, 1)
}

# upper and lower bounds of x
bounds = {
    "Fold": (-1, np.inf),
    "Transcritical": (-1, np.inf),
    "Pitchfork": (-1, 1)
}

# normal form equation in latex
eq = {
    "Fold": "f(x,y)=-y-x^2",
    "Transcritical": "f(x,y)=y x-x^2",
    "Pitchfork": "f(x)=y x+x^3"
}


st.header("Sample Path")

# bifurcation type
type = st.selectbox(
    "Bifurcation",
    options=types + ["Custom"]
)

if type == "Custom":
    type = st.text_input("Function Name")
    func = st.text_input("f(x, y)")

    # bounds
    lower = st.number_input("Lower Bound", max_value=0., value=-10., step=0.1)
    upper = st.number_input("Upper Bound", min_value=0., value=10., step=0.1)
    bounds[type] = (lower, upper)

    f[type] = lambda x, y: eval(func)
    x0[type] = st.text_input("x0")
    domain[type] = lambda y: bounds[type]
    xlim[type] = bounds[type]

    types.append(type)

# display normal form equation in latex
if type in eq:
    st.latex(eq[type])

# button to generate sample path
if st.button("Generate sample path"):
    path = sample_path(f[type], x0[type], y0, T, h, noise, scale, bounds=bounds[type])
    
    fig, ax = plt.subplots()
    ax.plot(path)
    ax.set_xlabel("y")
    ax.set_ylabel("x")
    st.pyplot(fig)


st.header("Sample Path Ensemble")

selected_types = st.multiselect("Bifurcations", options=types)

# number of sample paths to generate
runs = st.number_input("Number of Sample Paths Generated", min_value=1, value=30)

# number of sample paths to display
runs_displayed = st.slider("Number of Sample Paths Displayed", min_value=0, max_value=runs, value=min(30, runs))


if len(selected_types):
    st.subheader("Sample Path Trajectories")
    paths = {}

    for type in selected_types:
        paths[type] = ensemble_paths_cache(f[type], x0[type], y0, T, h, noise, scale, bounds=bounds[type], runs=runs)

    def plot_sample_paths(ax, type):
        plot_ensemble(
            ax,
            paths[type].iloc[:,:runs_displayed],
            equilibria=equilibria[type] if type in equilibria else []
        )
    
    st.pyplot(plot_grid(selected_types, "y", "x", plot_sample_paths))


    st.subheader("Percentage Escaped Trajectories")
    percent_escaped = {}

    for type in selected_types:
        percent_escaped[type] = escaped_trajectories_cache(paths[type], domain[type])

    def plot_escaped(ax, type):
        ax.plot(percent_escaped[type])
        ax.set_ylim((0, 100))

    st.pyplot(plot_grid(selected_types, "y", "% escaped", plot_escaped))


    st.subheader("Variance")

    estimated_var = {}
    analytical_var = {}

    for type in selected_types:
        estimated_var[type] = estimate_variance_cache(paths[type], domain[type])
        if type in pdf:
            analytical_var[type] = analytical_variance_cache(paths[type], domain[type], pdf[type])

    def plot_variance(ax, type):
        ax.plot(estimated_var[type])
        if type in analytical_var:
            ax.plot(analytical_var[type])

    st.pyplot(plot_grid(selected_types, "y", "Variance", plot_variance))


    st.subheader("Probability Distribution")

    # y value at time t
    yt = st.slider("yt", min_value=y0, max_value=min(0., yf))

    xt = {}  # distribution of x values at time t
    px = {}  # p(x) at given y

    for type in selected_types:
        xt[type] = distribution_cache(paths[type], yt)
        if type in pdf:
            px[type] = probability_density_cache(pdf[type], domain[type], yt)

    def plot_distribution(ax, type):
        ax.hist(xt[type], bins=bound(len(xt[type]) // 20, 5, 30), density=True)
        if type in px:
            ax.plot(px[type])
        ax.set_xlim(xlim[type])
    
    st.pyplot(plot_grid(selected_types, f"x(y={yt})", "Probability density", plot_distribution))


    st.subheader("Surface Plot")

    pxy = {}  # p(x,y)

    for type in selected_types:
        if type in pdf:
            pxy[type] = probability_grid(xlim[type], (y0, 0), pdf[type], domain[type])

    def plot_surface(ax, type):
        ax.contourf(*pxy[type], levels=np.arange(0, 15, 1), cmap=cm.coolwarm)
        plot_ensemble(ax, paths[type][paths[type].index < 0].iloc[:,:runs_displayed])

    st.pyplot(plot_grid([type for type in selected_types if type in pdf], "y", "x", plot_surface))