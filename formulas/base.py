from __future__ import print_function

import numpy as np

from sympy.utilities.lambdify import lambdify
from sympy import symbols, latex, Float, diff, \
    sympify, Mul, Add, Rational, Integer
from sympy.core.numbers import Pi

import matplotlib
from matplotlib.backends.backend_agg import \
    FigureCanvasAgg

import base64
from io import BytesIO
from collections import OrderedDict

from formulas.utils import \
    in_units, \
    latexify_units, \
    get_units, \
    check_type

NumberTypes = [int, float, Float, Rational, Integer, Pi]
NumberTypes = tuple(NumberTypes)
SymPyTypes = (Mul, Add)
ExpTypes = (Rational, Float, Integer)


class Formula(object):
    def __init__(self, formula, vars, params):
        self.function = None
        self.params = OrderedDict()
        self.param_values = OrderedDict()
        self.var_symbols = OrderedDict()
        for v in vars:
            if isinstance(v, str):
                v = symbols(v)
            self.var_symbols[str(v)] = v
        self.ndim = len(self.var_symbols.keys())
        for p in params:
            if isinstance(p, str):
                p = symbols(p)
            self.params[str(p)] = p
            self.param_values[str(p)] = None
        if isinstance(formula, str):
            formula = sympify(formula)
        self.formula = formula
        self.num_params = len(self.params)
        if self.num_params == 0:
            self.set_param_values()

    def __repr__(self):
        return self.formula.__repr__()

    def show(self): # pragma: no cover
        """
        Provide a pretty representation of the formula
        from within IPython.
        """
        from IPython.display import Latex
        return Latex(self.latex_representation(mode="inline"))

    def latex_representation(self, **kwargs):
        """
        Return the latex representation of the formula.
        """
        return latex(self.formula, **kwargs)

    def show_params(self): # pragma: no cover
        """
        Print out the parameters and their current values.
        """
        for k,v in self.param_values.items():
            print("%s = %s" % (k,v))

    def clear_param_values(self):
        """
        Set all of the parameter values to None.
        """
        for k in self.param_values:
            self.param_values[k] = None

    def set_param_values(self, **kwargs):
        """
        Set the values of one or more parameters.

        Examples
        --------
        >>> import yt.units as u
        >>> r_c = 50.*u.kpc
        >>> beta = 2./3.
        >>> rho_c = 1.0e8*u.Msun/u.kpc**3
        >>> density_profile = beta_model_profile()
        >>> density_profile.set_param_values(r_c=r_c, beta=beta, rho_c=rho_c)
        """
        for k,v in kwargs.items():
            if k in self.param_values:
                self.param_values[k] = v
            else:
                raise KeyError("Parameter %s is not in this formula!" % k)

    def copy(self):
        vars = list(self.var_symbols.values())
        params = list(self.params.values())
        if self.ndim == 1:
            f = Formula1D(self.formula, vars[0], params)
        elif self.ndim == 2:
            f = Formula2D(self.formula, vars[0], vars[1], params)
        else:
            f = Formula(self.formula, vars, params)
        f.set_param_values(**self.param_values)
        return f

    def __call__(self, **kwargs):
        if None in self.param_values.values():
            raise ValueError("Not all of the parameters are set!")
        if len(kwargs) != self.ndim:
            raise RuntimeError("Incorrect number of arguments provided! Must be %d!" % self.ndim)
        vars = [kwargs[var] for var in self.var_symbols]
        args = vars+list(self.param_values.values())
        if self.function is None:
            fargs = list(self.var_symbols.values())+list(self.params.values())
            self.function = lambdify(fargs, self.formula, modules="numpy")
        return self.function(*args)

    def _formula_op(self, op, other):
        param_values = self.param_values.copy()
        params = list(self.params.values())
        vars = list(self.var_symbols.values())
        if isinstance(other, Formula):
            formula = op(other.formula)
            vars += [sym for name, sym in other.var_symbols.items() if name not in self.var_symbols]
            for name, value in other.param_values.items():
                if name in param_values:
                    if value is not None and value != param_values[name]:
                        raise RuntimeError("The %s parameter is in both formulas, but has "
                                           "different values! Can't decide which one to "
                                           "use--choose one and set the other to None." % name)
                else:
                    params.append(other.params[name])
                    param_values[name] = value
        elif op is None and isinstance(other, SymPyTypes):
            formula = other
        elif isinstance(other, NumberTypes):
            formula = op(other)
        else:
            raise RuntimeError("Undefined operation between Formula and %s." % type(other))
        # The next couple of loops check for variables and parameters that have
        # been canceled out by subtract, division, or differentiation
        for var in list(vars):
            if var not in formula.free_symbols:
                vars.remove(var)
        for p in list(params):
            if p not in formula.free_symbols:
                params.remove(p)
                param_values.pop(str(p))
        ndim = len(vars)
        if ndim == 1:
            p = Formula1D(formula, vars[0], params)
        elif ndim == 2:
            p = Formula2D(formula, vars[0], vars[1], params)
        else:
            p = Formula(formula, vars, params)
        p.set_param_values(**param_values)
        return p

    def __add__(self, other):
        return self._formula_op(self.formula.__add__, other)

    def __radd__(self, other):
        return self._formula_op(self.formula.__radd__, other)

    def __sub__(self, other):
        return self._formula_op(self.formula.__sub__, other)

    def __rsub__(self, other):
        return self._formula_op(self.formula.__rsub__, other)

    def __neg__(self):
        formula = self.formula.__neg__()
        return self._formula_op(None, formula)

    def __pos__(self):
        formula = self.formula.__pos__()
        return self._formula_op(None, formula)

    def __mul__(self, other):
        return self._formula_op(self.formula.__mul__, other)

    def __rmul__(self, other):
        return self._formula_op(self.formula.__rmul__, other)

    def __truediv__(self, other):
        return self._formula_op(self.formula.__truediv__, other)

    def __pow__(self, other):
        if not isinstance(other, ExpTypes):
            if isinstance(other, int):
                a = Integer(other)
            else:
                a = Float(other)
        else:
            a = other
        return self._formula_op(self.formula.__pow__, a)

    __div__ = __truediv__

    def diff(self, var):
        r"""
        Return a new formula that is the derivative of the formula
        with respect to the given variable.
        """
        formula = diff(self.formula, self.var_symbols[var])
        return self._formula_op(None, formula)

    _unitless = None

    @property
    def unitless(self):
        if self._unitless is None:
            vars = list(self.var_symbols.values())
            params = list(self.params.values())
            if self.ndim == 1:
                uf = Formula1D(self.formula, vars[0], params)
            elif self.ndim == 2:
                uf = Formula2D(self.formula, vars[0], vars[1], params)
            else:
                uf = Formula(self.formula, vars, params)
            pvalues = {}
            for k,v in self.param_values.items():
                if hasattr(v,"units") or hasattr(v,"unit"):
                    pvalues[k] = float(v.value)
                else:
                    pvalues[k] = v
            uf.set_param_values(**pvalues)
            self._unitless = uf
        return self._unitless


class Formula1D(Formula):
    def __init__(self, formula, x, params):
        super(Formula1D, self).__init__(formula, [x], params)
        self.x = list(self.var_symbols.values())[0]

    def quick_plot(self, x_min, x_max, x_scale="linear", y_scale="linear",
                   res=200, filename=None, function_name="f", units=None): # pragma: no cover
        """
        Plot the formula.

        Parameters
        ----------
        x_min : unitful scalar quantity
            The mininum value of the "x" variable to plot.
        x_max : unitful scalar quantity
            The maximum value of the "x" variable to plot.
        x_scale : string
            The scaling for the x axis. Can be "linear" or "log".
        y_scale : string
            The scaling for the y axis. Can be "linear" or "log".
        res : integer
            The number of points to use in the plot. Default is 200.
        filename : str
            If set, save the plot to this filename.
        function_name : str
            The name of the function for the y-axis of the plot.
        units : str
            The units to convert the y-axis values to.

        Examples
        --------
        >>> import astropy.units as u
        >>> r_min = 0.1*u.kpc
        >>> r_max = 1000.*u.kpc
        >>> density_profile.quick_plot(r_min, r_max, x_scale="log")
        """
        from IPython.display import HTML, display
        matplotlib.rc("font", size=16, family="serif")
        fig = matplotlib.figure.Figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        arr = check_type(x_min)
        x = arr(np.linspace(x_min.value, x_max.value, num=res), get_units(x_min))
        y = arr(self(x))
        if units is not None:
            y = in_units(y, units)
        x_units = latexify_units(x)
        y_units = latexify_units(y)
        ax.plot(np.array(x), np.array(y))
        ax.set_xlabel(r"$\mathrm{%s}$ (" % self.x + x_units + ")")
        ax.set_ylabel(r"$\mathrm{%s(%s)}$ (" % (function_name, self.x) + y_units + ")")
        ax.set_xscale(x_scale)
        ax.set_yscale(y_scale)
        fig.tight_layout()
        if filename is not None:
            fig.savefig(filename)
        canvas = FigureCanvasAgg(fig)
        f = BytesIO()
        canvas.print_figure(f)
        f.seek(0)
        img = base64.b64encode(f.read()).decode()
        ret = r'<img style="max-width:100%%;max-height:100%%;" ' \
              r'src="data:image/png;base64,{0}"><br>'.format(img)
        display(HTML(ret))

    def __call__(self, x):
        kwargs = {str(self.x):x}
        return super(Formula1D, self).__call__(**kwargs)

class Formula2D(Formula):
    def __init__(self, formula, x, y, params):
        super(Formula2D, self).__init__(formula, [x,y], params)
        self.x, self.y = list(self.var_symbols.values())

    def quick_plot(self, x_min, x_max, y_min, y_max, x_scale="linear",
                   y_scale="linear", z_scale="linear", res_x=200,
                   res_y=200, filename=None, cmap=None, function_name="f",
                   units=None): # pragma: no cover
        """
        Plot the formula.

        Parameters
        ----------
        x_min : unitful scalar quantity
            The mininum value of the "x" variable to plot.
        x_max : unitful scalar quantity
            The maximum value of the "x" variable to plot.
        y_min : unitful scalar quantity
            The mininum value of the "y" variable to plot.
        y_max : unitful scalar quantity
            The maximum value of the "y" variable to plot.
        x_scale : string
            The scaling for the x axis. Can be "linear" or "log".
        y_scale : string
            The scaling for the y axis. Can be "linear" or "log".
        z_scale : string
            The scaling for the z axis. Can be "linear" or "log".
        res_x : integer
            The number of points to use in the plot along the x-axis.
            Default is 200.
        res_y : integer
            The number of points to use in the plot along the y-axis.
            Default is 200.
        filename : str
            If set, save the plot to this filename.
        cmap : str
            The colormap to use when making the plot.
        function_name : str
            The name of the function for the colormap.
        units : str
            The units to convert the colormap values to.

        Examples
        --------
        >>> import yt.units as u
        >>> g_x = gaussian(x="v_x", A="A_x", mu="mu_x", sigma="sigma_x")
        >>> g_y = gaussian(x="v_y", A="A_y", mu="mu_y", sigma="sigma_y")
        >>> g = g_x*g_y
        >>> g.set_param_values(A_x=1.0, A_y=1.0, mu_x=0*u.km/u.s, mu_y=0*u.km/u.s,
        ...                    sigma_x=200*u.km/u.s, sigma_y=100*u.km/u.s)
        >>> g.quick_plot(-300*u.km/u.s, 300*u.km/u.s, -300*u.km/u.s, 300*u.km/u.s,)
        """
        from IPython.display import display, HTML
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        matplotlib.rc("font", size=16, family="serif")
        fig = matplotlib.figure.Figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        arr = check_type(x_min)
        x = np.linspace(x_min.value, x_max.value, num=res_x)
        y = np.linspace(y_min.value, y_max.value, num=res_y)
        xx, yy = np.meshgrid(x, y)
        xx = arr(xx, get_units(x_min))
        yy = arr(yy, get_units(y_min))
        vars = {str(self.x):xx,str(self.y):yy}
        z = arr(self(**vars))
        if units is not None:
            z = in_units(z, units)
        x_units = latexify_units(xx)
        y_units = latexify_units(yy)
        z_units = latexify_units(z)
        extent = (x[0], x[-1], y[0], y[-1])
        im = ax.imshow(z, extent=extent, cmap=cmap)
        ax.set_xlabel(r"$\mathrm{%s}$ (" % self.x + x_units + ")")
        ax.set_ylabel(r"$\mathrm{%s}$ (" % self.y + y_units + ")")
        ax.set_xscale(x_scale)
        ax.set_yscale(y_scale)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label(r"$\mathrm{%s(%s,%s)}$ (" % (function_name, self.x, self.y) + z_units + ")")
        cbar.ax.set_yscale(z_scale)
        fig.tight_layout()
        if filename is not None:
            fig.savefig(filename)
        canvas = FigureCanvasAgg(fig)
        f = BytesIO()
        canvas.print_figure(f)
        f.seek(0)
        img = base64.b64encode(f.read()).decode()
        ret = r'<img style="max-width:100%%;max-height:100%%;" ' \
              r'src="data:image/png;base64,{0}"><br>'.format(img)
        display(HTML(ret))

def variable(x):
    x = symbols(x)
    return Formula1D(x, x, [])
