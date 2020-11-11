"""The Wavefunction class and its attributes."""

import numpy as np
import pyfftw
import numexpr as ne
import pytalises.propagator


class Wavefunction:
    """
    Class describing wave function.

    Parameters
    ------------------
    psi : string list of strings
        Each string describes the initial amplitude in r-space of
        the wave function. The number of elements is the number
        of internal degrees of freedom.
        Additional variables that are used in psi can
        be defined in the optional parameter variables.
        Predefined variables are the spatial coordinates x,y,z
        and time t.
    number_of_grid_points : tuple of ints
        Tuple that defines the number of grid points (nX,nY,nZ)
        of the wave function
    spatial_ext : list of tuples
        The supplied values define the boundary positions of the grid
        and thus define the actual coordinate system.
    t0: float, optional
        Internal time of wave function. Default is 0.0.
    m : float, optional
        Mass of particle described by the wavefunction.
        Default is 1.054571817e-34 (numerically equal to hbar).
    variables : dict
        Dictionary of additionaly used variables in the definition of
        the wave function in psi.
        Predefined variables are the spatial coordinates x,y,z
        and time t.
    normalize_const : float, optional
        Normalizes the wave function such that the integral of |Psi|^2
        over all internal and external degrees of freedom equals
        `normalize_const`

    Attributes
    ------------------
    num_int_dim : int
        The number of internal degrees of freedom
    num_ex_dim : int
        The number of external degrees of freedom
    r : list of 1d arrays
        The arrays are the evenly spaced spatial coordinates as defined
        through definition of spatial_ext and number_of_grid_points
    k : list of 1d arrays
        The arrays are the evenly spaced inerse spatial coordinates as defined
        through definition of spatial_ext and number_of_grid_points

    Examples
    ------------------
    Wavefunction with two interal states
    where the first state is gaussian
    distributed in 1d r-space and the second
    state is not occupied at all.

    >>> from pytalises import Wavefunction
    >>> psi = Wavefunction(["exp(-((x-x0)/a0)**2)", "0.0"],
        (16,), [(-2,2),], variables={'a0':1/2, 'x0':0})
    >>> print(psi.num_int_dim)
    2
    >>> print(psi.num_ext_dim)
    1
    >>> print(psi.amp)
    [[1.12535175e-07+0.j 0.00000000e+00+0.j]
    [6.03594712e-06+0.j 0.00000000e+00+0.j]
    [1.83289361e-04+0.j 0.00000000e+00+0.j]
    [3.15111160e-03+0.j 0.00000000e+00+0.j]
    [3.06707930e-02+0.j 0.00000000e+00+0.j]
    [1.69013315e-01+0.j 0.00000000e+00+0.j]
    [5.27292424e-01+0.j 0.00000000e+00+0.j]
    [9.31358402e-01+0.j 0.00000000e+00+0.j]
    [9.31358402e-01+0.j 0.00000000e+00+0.j]
    [5.27292424e-01+0.j 0.00000000e+00+0.j]
    [1.69013315e-01+0.j 0.00000000e+00+0.j]
    [3.06707930e-02+0.j 0.00000000e+00+0.j]
    [3.15111160e-03+0.j 0.00000000e+00+0.j]
    [1.83289361e-04+0.j 0.00000000e+00+0.j]
    [6.03594712e-06+0.j 0.00000000e+00+0.j]
    [1.12535175e-07+0.j 0.00000000e+00+0.j]]
    >>> print(psi.r)
    [array([-2.        , -1.73333333, -1.46666667, -1.2       , -0.93333333,
       -0.66666667, -0.4       , -0.13333333,  0.13333333,  0.4       ,
        0.66666667,  0.93333333,  1.2       ,  1.46666667,  1.73333333,
        2.        ]), array([0.]), array([0.])]
    """

    def __init__(
        self,
        psi,
        number_of_grid_points,
        spatial_ext,
        t0=0.0,
        m=1.054571817e-34,
        variables={},
        normalize_const=None,
    ):
        """Initialize Wavefunction."""
        if type(psi) is not list and type(psi) is str:
            psi = [psi]
        if type(spatial_ext) is not list and type(spatial_ext) is tuple:
            spatial_ext = [spatial_ext]
        self.num_int_dim = len(psi)
        self.num_ext_dim = sum([1 for n in number_of_grid_points if n > 0])
        assert isinstance(number_of_grid_points, tuple)
        for _ in range(3 - len(number_of_grid_points)):
            number_of_grid_points += (1,)
        self.nX, self.nY, self.nZ = number_of_grid_points
        self.number_of_grid_points = number_of_grid_points
        self.spatial_ext = spatial_ext
        assert self.num_ext_dim == len(spatial_ext)
        for _ in range(3 - len(spatial_ext)):
            spatial_ext += ((0, 0),)
        self.axes = tuple(np.where(np.array(number_of_grid_points) > 1, 1, 0))
        r = []
        Delta_r = []
        k = []
        delta_k = []
        for i, spatial_ext_tuple in enumerate(spatial_ext):
            r_min = spatial_ext_tuple[0]
            r_max = spatial_ext_tuple[1]
            r.append(np.linspace(r_min, r_max, num=number_of_grid_points[i]))
            if self.axes[i] == 0:
                Delta_r.append(np.nan)
                delta_k.append(np.nan)
                k.append(0.0)
            else:
                Delta_r.append(r_max - r_min)
                delta_k.append(2 * np.pi / (r_max - r_min))
                k.append(
                    np.fft.fftfreq(number_of_grid_points[i])
                    * 2
                    * np.pi
                    * number_of_grid_points[i]
                    / Delta_r[i]
                )

        self._r = r
        self.Delta_r = Delta_r
        self.delta_r = [
            Delta / self.number_of_grid_points[i]
            for i, Delta in enumerate(self.Delta_r)
        ]
        self.delta_k = delta_k
        self._k = k
        self.rmesh = np.meshgrid(*r, indexing="ij")
        self.kmesh = np.meshgrid(*k, indexing="ij")
        self._amp = pyfftw.empty_aligned(
            number_of_grid_points + (self.num_int_dim,), dtype="complex128", order="C"
        )
        self.psi = psi
        self.t = t0
        self.m = m
        self.alpha = 1.054571817e-34 / (2 * self.m)
        self.default_var_dict = {
            "alpha": self.alpha,
            "x": self.rmesh[0],
            "y": self.rmesh[1],
            "z": self.rmesh[2],
        }
        for i in range(self.num_int_dim):
            self.default_var_dict["psi" + str(i)] = self._amp[:, :, :, i]
        self.variables = variables
        for i in range(self.num_int_dim):
            self._amp[:, :, :, i] = ne.evaluate(
                self.psi[i],
                local_dict={**self.default_var_dict, **self.variables},
                order="C",
            )
        self.normalize_const = normalize_const
        if normalize_const is not None:
            self.normalize_to(normalize_const)
        self.construct_FFT()

    def construct_FFT(
        self,
        num_of_threads=1,
        FFTWflags=(
            "FFTW_ESTIMATE",
            "FFTW_DESTROY_INPUT",
        ),
    ):
        """Construct pyfftw bindings."""
        axes = tuple(i for i in range(self.num_ext_dim))
        pyfftw.config.NUM_THREADS = num_of_threads
        self.fft = pyfftw.FFTW(
            self._amp,
            self._amp,
            axes=axes,
            direction="FFTW_FORWARD",
            threads=num_of_threads,
            flags=FFTWflags,
        )
        self.ifft = pyfftw.FFTW(
            self._amp,
            self._amp,
            axes=axes,
            direction="FFTW_BACKWARD",
            threads=num_of_threads,
            flags=FFTWflags,
        )

    @property
    def r(self):
        """(list of) array of the wave function position axes."""
        if self.num_ext_dim == 1:
            return self._r[self.axes.index(1)]
        else:
            return self._r

    @property
    def k(self):
        """(list of) array of the wave function position axes."""
        if self.num_ext_dim == 1:
            return self._k[self.axes.index(1)]
        else:
            return self._k

    @property
    def amp(self):
        """Ndarray of the wave function amplitudes."""
        return np.squeeze(self._amp)

    def exp_pos(self, axis=None):
        """
        Calculate the expected position on given axis.

        Calculates the mean position of Psi on chosen axis.
        Axes 0,1,2 correspond to x,y,z. The other two axes
        are traced out. If no axis iv given returns array
        of mean position of all external degrees of freedom.
        """
        if axis is None:
            exp_pos = np.empty((self.num_ext_dim))
            for i in range(self.num_ext_dim):
                exp_pos[i] = self.exp_pos(i)
        else:
            axes_to_trace = [0, 1, 2]
            axis = axes_to_trace.pop(axis)
            psi_sq_amp = np.power(np.abs(self._amp), 2)
            traced_out_psi = np.sum(psi_sq_amp, axis=tuple(axes_to_trace))
            exp_pos = np.einsum("ri,r->", traced_out_psi, self._r[axis])
            exp_pos *= np.prod(self.delta_r, where=np.where(self.axes, True, False))
        return exp_pos

    def normalize_to(self, n_const):
        """
        Normalize the wave function.

        Normalizes the wave function such that the integral
        of |Psi|^2 over all internal and external states
        equals n_const
        """
        # Calculate |Psi|^2 over all internal and external states
        s = np.einsum("xyzi,xyzi->", self._amp, np.conjugate(self._amp))
        # Mulitply with product of infinitesimal volumes dx*dy*dz
        # while ignoring nonextisten dimensions
        s *= np.prod(self.delta_r, where=np.where(self.axes, True, False))
        self._amp *= np.sqrt(n_const / s)

    def state_occupation(self, nth_state=None):
        """
        Return occupation number of nth internal state.

        Evaluates the spatial integral over |Psi|^2 for
        the nth internal state. If none is given a vector
        of the occupation number of all internal states
        is returned.
        """
        if nth_state is None:
            state_occupation = np.empty((self.num_int_dim,))
            for i in range(self.num_int_dim):
                state_occupation[i] = self.state_occupation(i)
        else:
            state_occupation = np.sum(np.abs(self.amp[:, nth_state]) ** 2) * np.prod(
                self.delta_r, where=np.where(self.axes, True, False)
            )
        return state_occupation

    def freely_propagate(
        self,
        num_time_steps,
        delta_t,
        num_of_threads=1,
        FFTWflags=(
            "FFTW_ESTIMATE",
            "FFTW_DESTROY_INPUT",
        ),
    ):
        """
        Propagates the Wavefunction object in time with V=0.

        Function that can propagate the wavefunction if no potential
        is present.

        Parameters
        ------------------
        num_time_steps : int
            Number of times the wavefunction is propagated by time delta_t
            using the Split-Steo Fourier method.
        delta_t : float
            Time increment the wavefunction is propagated in one time step.
        num_of_threads : int, optional
            Number of threads uses for calculation. Default is 1.
        FFTWflags : tuple of strings
            Options for FFTW planning [1]. Default is
            ('FFTW_ESTIMATE', 'FFTW_DESTROY_INPUT',).

        References
        --------
        [1] http://www.fftw.org/fftw3_doc/Planner-Flags.html
        """
        pytalises.propagator.freely_propagate(
            self, num_time_steps, delta_t, num_of_threads, FFTWflags
        )

    def propagate(self, potential, num_time_steps, delta_t, **kwargs):
        """
        Propagates the Wavefunction object in time.

        Function that propagates the wavefunction using a
        Split-Step Fourier method [1].

        Parameters
        ------------------
        potential : string list of strings
            This list contains the matrix elements of the potential term V
            in string format. If the potential has nondiagonal elements
            (see optional parameter diag) earch elements represents
            one matrix element of the lower triangular part of V.
            For example a 3x3 potential with nondiagonal elements would be
            of form potential=[H00, H10, H20, H11, H21, H22].
            If the potential term is supposed to have only diagonal elements
            (diag=True), the potential parameter for a 3x3 potential would
            look like potential=[H00,H11,H22].
        num_time_steps : int
            Number of times the wavefunction is propagated by time delta_t
            using the Split-Steo Fourier method.
        delta_t : float
            Time increment the wavefunction is propagated in one time step.
        variables : dict, optional
            Dictionary containing values for variables you might have used
            in potential
        diag : bool , optional
            If true, no numerical diagonalization has to be invoked in order
            to calculate time-propagation as nondiagonal elements are ommited.
            This makes the computation much faster. Default is False.
        num_of_threads : int, optional
            Number of threads uses for calculation. Default behaviour
            is to use all threads available.
        FFTWflags : tuple of strings
            Options for FFTW planning [2]. Default is
            ('FFTW_ESTIMATE', 'FFTW_DESTROY_INPUT',).

        References
        --------
        [1] https://en.wikipedia.org/wiki/Split-step_method
        [2] http://www.fftw.org/fftw3_doc/Planner-Flags.html
        """
        pytalises.propagator.propagate(
            self, potential, num_time_steps, delta_t, **kwargs
        )
