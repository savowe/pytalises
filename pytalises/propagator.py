"""Module containing functions that help propagating the Wavefunction class."""
from numba import jit, prange, set_num_threads
from numpy.linalg import eigh
import numexpr as ne
import numpy as np
import pytalises.wavefunction


def propagate(psi, potential, num_time_steps, delta_t, **kwargs):
    """
    Propagates a Wavefunction object in time.

    Function that propagates the wavefunction using a
    Split-Step Fourier method [1].

    Parameters
    ------------------
    psi : Wavefunction
        The Wavefunction object the Propagator class acts on
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
    U = Propagator(psi, potential, **kwargs)
    U.kinetic_prop(delta_t / 2)
    U.potential_prop(delta_t)
    for _ in range(num_time_steps - 1):
        U.kinetic_prop(delta_t)
        U.potential_prop(delta_t)
    U.kinetic_prop(delta_t / 2)


def freely_propagate(
    psi,
    num_time_steps,
    delta_t,
    num_of_threads=1,
    FFTWflags=(
        "FFTW_ESTIMATE",
        "FFTW_DESTROY_INPUT",
    ),
):
    """
    Propagates a Wavefunction object in time with V=0.

    Function that can propagate the wavefunction if no potential
    is present.

    Parameters
    ------------------
    psi : Wavefunction
        The Wavefunction object the Propagator class acts on
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
    U = Propagator(
        psi,
        potential=["0"] * psi.num_int_dim,
        diag=True,
        num_of_threads=num_of_threads,
        FFTWflags=FFTWflags,
    )
    for _ in range(num_time_steps):
        U.kinetic_prop(delta_t)


class Propagator:
    """
    Class for propagating instances of the Wavefunction class.

    Parameters
    ------------------
    psi : Wavefunction
        The Wavefunction object the Propagator class acts on
    potential : list of strings
        This list contains the matrix elements of the potential term V
        in string format. If the potential has nondiagonal elements
        (see optional parameter diag) earch elements represents
        one matrix element of the lower triangular part of V.
        For example a 3x3 potential with nondiagonal elements would be
        of form potential=[H00, H10, H20, H11, H21, H22].
        If the potential term is supposed to have only diagonal elements
        (diag=True), the potential argument for a 3x3 potential would
        look like potential=[H00,H11,H22].
    variables : dict, optional
        Dictionary containing values for variables you might have used
        in potential
    diag : bool , optional
        If true, no numerical diagonalization has to be invoked in order
        to calculate time-propagation. Default is False.
    num_of_threads : int, optional
        Number of threads uses for calculation. Default is 1.
    FFTWflags : tuple of strings
        Options for FFTW planning [1]. Default is
        ('FFTW_ESTIMATE', 'FFTW_DESTROY_INPUT',).

    References
    --------
    [1] http://www.fftw.org/fftw3_doc/Planner-Flags.html
    """

    def __init__(
        self,
        psi,
        potential,
        variables={},
        diag=False,
        num_of_threads=1,
        FFTWflags=(
            "FFTW_ESTIMATE",
            "FFTW_DESTROY_INPUT",
        ),
    ):
        """Initialize the propagator."""
        self.psi = psi
        self.v = self.Potential(potential, variables, diag)
        assert isinstance(psi, pytalises.wavefunction.Wavefunction)
        assert self.v.num_int_dim == self.psi.num_int_dim
        assert self.psi._amp.shape[-1] == self.psi.num_int_dim
        self.V_eval_array = np.zeros(
            psi.number_of_grid_points + (psi.num_int_dim, psi.num_int_dim),
            order="C",
            dtype="complex128",
        )
        self.V_eval_eigval_array = np.zeros(
            psi.number_of_grid_points + (psi.num_int_dim,),
            order="C",
            dtype="complex128",
        )
        self.num_of_threads = num_of_threads
        set_num_threads(num_of_threads)
        ne.set_num_threads(num_of_threads)
        self.psi.construct_FFT(num_of_threads, FFTWflags)
        if self.v.diag is True:
            self.prop_method = self.diag_potential_prop
        else:
            self.prop_method = self.nondiag_potential_prop
        # Check if potential is static and if that is the case
        # precompute the potential grid V(x,y,z) to use it
        # for all following calculations
        if self.v.static is True:
            self.eval_V()
            if self.v.diag is False:
                get_eig(self.V_eval_array, self.V_eval_eigval_array)

    def potential_prop(self, delta_t):
        """
        Wrap function that calculates exp(i*V(x,y,z)/hbar*delta_t)*Psi(x,y,z).

        This can be either nondiag_potential_prop or diag_potential_prop.
        """
        self.prop_method(delta_t)

    def nondiag_potential_prop(self, delta_t):
        """
        Calculate exp(i*V/hbar*delta_t)*Psi using numerical diagonalization.

        This method has to be used if the potential mmatrix has nondiagonal
        elements.
        """
        if self.v.static is False:
            print("potential is time-dependent")
            self.eval_V()
            get_eig(self.V_eval_array, self.V_eval_eigval_array)
        np.einsum(
            "xyzij,xyzj,xyzkj,xyzk->xyzi",
            self.V_eval_array,
            ne.evaluate(
                "exp(-1j*eigval*delta_t)",
                local_dict={"eigval": self.V_eval_eigval_array, "delta_t": delta_t},
            ),
            np.conjugate(self.V_eval_array),
            self.psi._amp,
            out=self.psi._amp,
            optimize="optimal",
            order="C",
        )

    def diag_potential_prop(self, delta_t):
        """
        Calculate exp(i*V/hbar*delta_t)*Psi by simple matrix multiplication.

        This method is used if the potential matrix V is diagonal. This is
        much faster than `nondiag_potential_prop` and should be used if
        possible.
        """
        if self.v.static is False:
            self.eval_diag_V()
        np.einsum(
            "xyzii,xyzi->xyzi",
            ne.evaluate(
                "exp(-1j*V*delta_t)",
                local_dict={"V": self.V_eval_array, "delta_t": delta_t},
            ),
            self.psi._amp,
            out=self.psi._amp,
            optimize="optimal",
            order="C",
        )

    def kinetic_prop(self, delta_t):
        """
        Perform time propagation in k-space.

        Transforms the Wavefunction into k-space,
        calculates exp(i*hbar/(2m)*k**2*delta_t)*Psi(kx,ky,kz)
        and transforms it back into r-space.
        """
        self.psi.fft()
        np.einsum(
            "xyz,xyzi->xyzi",
            ne.evaluate(
                "exp(-1j*alpha*delta_t*(kx**2+ky**2+kz**2))",
                local_dict={
                    "kx": self.psi.kmesh[0],
                    "ky": self.psi.kmesh[1],
                    "kz": self.psi.kmesh[2],
                    "alpha": self.psi.alpha,
                    "delta_t": delta_t,
                },
                order="C",
            ),
            self.psi._amp,
            out=self.psi._amp,
            optimize="optimal",
        )
        self.psi.ifft()
        self.psi.t += delta_t

    def eval_V(self):
        """
        Evalutes V on the whole spatial grid.

        The result is saved in Propagator.V_eval_array.
        """
        k = 0
        for i in range(self.psi.num_int_dim):
            for j in range(i, self.psi.num_int_dim):
                self.V_eval_array[:, :, :, j, i] = ne.evaluate(
                    self.v.potential_strings[k],
                    local_dict={**self.v.variables, **self.psi.default_var_dict},
                    global_dict={"t": self.psi.t},
                    order="C",
                )
                k += 1

    def eval_diag_V(self):
        """
        Evalutes diagonal elements of V on the whole spatial grid.

        The result is saved in Propagator.V_eval_array.
        """
        for i in range(self.psi.num_int_dim):
            self.V_eval_array[:, :, :, i, i] = ne.evaluate(
                self.v.potential_strings[i],
                local_dict={**self.v.variables, **self.psi.default_var_dict},
                global_dict={"t": self.psi.t},
                order="C",
            )

    class Potential:
        """Simple class for collecting information about the potential."""

        def __init__(self, potential_string, variables={}, diag=False):
            """Initialize Potential."""
            if type(potential_string) is not list and type(potential_string) is str:
                self.potential_strings = [potential_string]
            else:
                self.potential_strings = potential_string
            self.num_v = len(self.potential_strings)
            self.variables = variables
            # Check if potential is static in time
            for pot_string in self.potential_strings:
                potential_nex = ne.NumExpr(pot_string)
                try:
                    potential_nex.input_names.index("t")
                    self.static = False
                except ValueError:
                    self.static = True
                if self.static is False:
                    break
            # Check if potential is linear (independent of psi).
            # If it depends on psi, it also depends on t and
            # is therefore not static
            for i in range(len(self.potential_strings)):
                for pot_string in self.potential_strings:
                    potential_nex = ne.NumExpr(pot_string)
                    try:
                        potential_nex.input_names.index("psi" + str(i))
                        self.linear = False
                    except ValueError:
                        self.linear = True
                    if self.linear is False:
                        self.static = False
                        break
                if self.linear is False:
                    self.static = False
                    break
            # Check if number of matrix elements given matches
            # the number diagonal or nondiagonal hermitian matrix
            # elements. In the case of a diagonal matrix the number
            # of matrix elements is equal to the number of internal
            # states. If it is nondiagonal one gives the lower
            # triangular part of V.
            self.diag = diag
            if diag is False:
                self.num_int_dim = 1 / 2 * (np.sqrt(8 * self.num_v + 1) - 1)
                assert (
                    self.num_int_dim.is_integer()
                ), "Number of potential matrix elements incorrect"
                self.num_int_dim = int(self.num_int_dim)
            if diag is True:
                self.num_int_dim = len(self.potential_strings)


@jit(nopython=True, parallel=True, nogil=True, fastmath=True)
def get_eig(matrices, eigvals):
    """
    Calculate eigenvectors and eigenvalues of matrices in array.

    JIT-compiled function that calculates the eigenvectors and
    eigenvalues of input array M in parallel using numba.
    The resulting eigenvectors are stored in the input matrix
    and the eigenvalues in the array eigvals.

    Parameters
    ------------------
    M : 3d array of (NxN) arrays
    eigvals : 3d array of 1d arrays with N elements
    """
    nX, nY, nZ = matrices.shape[:3]
    for i in prange(nX):
        for j in prange(nY):
            for k in prange(nZ):
                eigvals[i, j, k, :], matrices[i, j, k, :, :] = eigh(
                    matrices[i, j, k, :, :]
                )
