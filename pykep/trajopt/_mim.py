import pykep as _pk
import heyoka as _hy
import numpy as _np
from scipy.optimize import root as _root


def _mim_shooting(x, ta, posvel0, posvelf, tof):
    """Terminal residuals of the MIM shooting problem, x = [lx, ly, lz, lvx, lvy, lvz, lm, m]."""
    ta.time = 0.0
    ta.state[:6] = posvel0
    ta.state[6] = x[7]
    ta.state[7:14] = x[:7]
    ta.propagate_until(tof)
    # Rendezvous, lm(tf) = 0 (squares the system) and |lambdas|^2 = 1
    return _np.r_[ta.state[:6] - posvelf, ta.state[13], x[:7] @ x[:7] - 1.0]


def _mim_jacobian(x, ta_var, ic_var, posvel0, tof):
    """Jacobian of :func:`_mim_shooting` from the variational equations w.r.t. [lx, ly, lz, lvx, lvy, lvz, lm, m]."""
    ta_var.time = 0.0
    ta_var.state[:6] = posvel0
    ta_var.state[6] = x[7]
    ta_var.state[7:14] = x[:7]
    ta_var.state[14:] = ic_var
    ta_var.propagate_until(tof)
    sensitivities = [
        ta_var.state[ta_var.get_vslice(order=1, component=i)] for i in (0, 1, 2, 3, 4, 5, 13)
    ]
    return _np.vstack(sensitivities + [_np.r_[2.0 * x[:7], 0.0]])


def mim_from_hop(pl_s, pl_t, when_s, when_f, Tmax, veff, MASS=1000.0, max_trials=10):
    """mim_from_hop(pl_s, pl_t, when_s, when_t, Tmax, veff, MASS=1000.0, max_trials=10)

    The Maximum Initial Mass from a hop transfer. A hop transfer is fixed time randevouz between two planets.
    This function is slower then its approximations :class:`~pykep.mima_from_hop` and :class:`~pykep.mima2_from_hop`,
    if succesfull, returns the exact solution to the problem. Behind the scene Potryagin Maximum Principle
    is applied to the problem and the TPBVP is solved by single shooting using :func:`scipy.optimize.root`
    (MINPACK hybrid method) and heyoka.

    Args:

        *pl_s* (:class:`~pykep.planet`): the source planet

        *pl_t* (:class:`~pykep.planet`): the target planet

        *when_s* (:class:`~pykep.epoch`): the start epoch

        *when_f* (:class:`~pykep.epoch`): the final epoch

        *Tmax* (:class:`float`, optional): Maximum spacecraft thrust.

        *veff* (:class:`float`, optional): Isp*G0.

        *MASS* (:class:`float`, optional): Units to use to scale the mass intrnally. It also influences the range of the random
        initial guesses for the MIM.           

        *max_trials* (:class:`int`, optional): Maximum number of random initial guesses tried to find the MIM. If the algorithm does not converge,
        it will return None.                        

    Returns:
        :class:`float`: The Maximum Initial Mass

    Examples:   

        >>> import pykep as pk
        >>> ... # assuming to have two planets pl_s and pl_t
        >>> when_s = pk.epoch(0.0)
        >>> when_f = pk.epoch(100.0)
        >>> Tmax = 0.6
        >>> veff = 3000*pk.G0
        >>> mima = pk.mim_from_hop(pl_s, pl_t, when_s, when_f, Tmax, veff)
        >>> print("Maximum initial mass:", mima, "kg")  
    """
    # Non-dimensional units
    L = pl_s.elements()[0]
    TIME = pl_s.period()
    VEL = L / TIME
    ACC = VEL / TIME
    mu = pl_s.get_mu_central_body() / (L**3 / TIME**2)
    c1 = Tmax / (MASS * ACC)
    c2 = veff / VEL
    tof = (when_f.mjd2000 - when_s.mjd2000) * _pk.DAY2SEC / TIME

    # Boundary conditions
    r0, v0 = pl_s.eph(when_s)
    rf, vf = pl_t.eph(when_f)
    posvel0 = _np.r_[r0, v0] / _np.r_[[L] * 3, [VEL] * 3]
    posvelf = _np.r_[rf, vf] / _np.r_[[L] * 3, [VEL] * 3]

    # Taylor integrators: the variational one has sensitivities w.r.t. the initial costates and mass
    ta = _pk.ta.get_pc(1e-16, _pk.optimality_type.TIME)
    dyn = _pk.ta.pc_dyn(_pk.optimality_type.TIME)
    args = _hy.make_vars("lx", "ly", "lz", "lvx", "lvy", "lvz", "lm", "m")
    ta_var = _hy.taylor_adaptive(
        _hy.var_ode_sys(sys=dyn, args=args, order=1), tol=1e-8, compact_mode=True
    )
    ic_var = ta_var.state[14:].copy()
    ta.pars[:3] = [mu, c1, c2]
    ta_var.pars[:3] = [mu, c1, c2]

    # Single shooting from random initial guesses
    rng = _np.random.default_rng()
    for i in range(max_trials):
        guess = _np.r_[rng.uniform(-1.0, 1.0, 6), rng.uniform(0.0, 1.0), rng.uniform(0.0, 20.0)]
        result = _root(
            lambda x: _mim_shooting(x, ta, posvel0, posvelf, tof),
            guess,
            jac=lambda x: _mim_jacobian(x, ta_var, ic_var, posvel0, tof),
            method="hybr",
            tol=1e-10,
            options={"factor": 1.0},
        )
        residual = _mim_shooting(result.x, ta, posvel0, posvelf, tof)
        # A root is accepted only if it is physical: positive initial and final mass
        if _np.linalg.norm(residual) < 1e-6 and result.x[7] > 0.0 and ta.state[6] > 0.0:
            return result.x[7] * MASS
    print(
        "Could not converge to find the MIM. tof might be too long or more trials may be needed (max_trials)"
    )
