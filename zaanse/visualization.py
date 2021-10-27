import numpy as np
import matplotlib.pyplot as plt


def plot_scalar_field(field, levels=10, title=""):
    """Plots a scalar field on a grid (phi, z).

    Our matrices save information such that z=0 corresponds to the upper row of the matrix.
    To create the plot, the matrix is mirrored vertically, i.e., the row at z=0 is below.

    Parameters
    ----------
    field : M by N array,
        The scalar field.
    levels : int, default=10
        The number of colors which will be used to draw the plot
    title : string, default=""
        The title of the plot
    """

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    height, width = field.shape
    x = np.linspace(0, 1, width, endpoint=False)
    y = np.arange(0, height)
    X, Y = np.meshgrid(x, y)

    contour_plot = ax.contourf(X, Y, field, levels, cmap=plt.cm.viridis)

    ax.set_xlabel(r"Angle $\phi$ [rad]")
    ax.set_ylabel(r"Height index")
    ax.set_title(title)

    # represent x-axis ticks in units of pi
    ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))

    fig.colorbar(contour_plot)

    plt.savefig(f"plots/{title}.png", dpi=300, bbox_inches="tight")


def plot_vector_field(field_x, field_y, scalar_field=None, density=(1, 2), title=""):
    """Plots a vector field on a grid (phi, z).

    Parameters
    ----------
    field_x : M by N array,
        The x component of the velocity.
    field_y : M by N array,
        The y component of the velocity.
    scalar_field : M by N array, optional
        Plots a scalar field in addition to the stream lines. For example,
        the norm of the vectors.
    density : float or tuple of two values
        Controls the closeness of streamlines. When density = 1, the domain
        is divided into a 30x30 grid. density linearly scales this grid.
        Each cell in the grid can have, at most, one traversing streamline.
        For different densities in each direction, use a tuple (density_x, density_y).
    title : string, default=""
        The title of the plot
    """
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))

    assert field_x.shape == field_y.shape
    height, width = field_x.shape
    x = np.linspace(0, 1, width, endpoint=False)
    y = np.arange(0, height)
    X, Y = np.meshgrid(x, y)

    ax.set_ylim([0, height - 1])

    linecolor = "black"
    # plot a scalar field in the background
    if type(scalar_field) != type(None):
        contour_plot = ax.contourf(
            X,
            Y,
            scalar_field,
            cmap=plt.cm.inferno,
            alpha=0.9,
            levels=30,
            linestyles=None,
        )
        fig.colorbar(contour_plot)
        linecolor = "white"

    ax.streamplot(
        X,
        Y,
        field_x,
        field_y,
        linewidth=0.5,
        density=density,
        color=linecolor,
        arrowsize=0.7,
    )

    ax.set_xlabel(r"Angle $\phi$ [rad]")
    ax.set_ylabel(r"Height index")
    ax.set_title(title)

    # represent x-axis ticks in units of pi
    ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))

    plt.savefig(f"plots/{title}.png", dpi=300, bbox_inches="tight")


# this solution is taken from https://stackoverflow.com/a/53586826

def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    def _multiple_formatter(x, pos):
        den = denominator
        num = int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
    return _multiple_formatter

class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex

    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)

    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))
