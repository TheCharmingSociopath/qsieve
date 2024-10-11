from magic_state_factory import MagicStateFactory
import mpmath
from mpmath import mp
mp.prec = 256
from scipy import optimize
from definitions import (
    z,
    one,
    projx,
    kron,
    trace,
    apply_rot,
    plog,
    ideal15to1,
    storage_x_4,
    storage_z_4,
    init4qubit,
    ideal8toCCZ,
)
from onelevel15to1 import one_level_15to1_state
from twolevel15to1 import two_level_15to1_state


def three_level_8toccz_state(
    pphys: float | mpmath.mpf,
    dx: int,
    dz: int,
    dm: int,
    dx2: int,
    dz2: int,
    dm2: int,
    dx3: int,
    dz3: int,
    dm3: int,
    nl1: int,
    nl2: int,
) -> (mpmath.matrix, float):
    """
    Generates the output-state density matrix of the (15-to-1)x(15-to-1)x(8-to-CCZ) protocol
    """

    pphys = mp.mpf(pphys)

    # Introduce shorthand notation for logical error rate with distances dx3/dz3/dm3
    px3 = plog(pphys, dx3)
    pz3 = plog(pphys, dz3)
    pm3 = plog(pphys, dm3)

    # Compute pl2, the output error of level-2 states
    out2, l1time = two_level_15to1_state(pphys, dx, dz, dm, dx2, dz2, dm2, nl1)
    pfail2 = 1 - trace(kron(one, projx, projx, projx, projx) * out2).real

    outpostsel2 = (
        (1 / (1 - pfail2))
        * kron(one, projx, projx, projx, projx)
        * out2
        * kron(one, projx, projx, projx, projx).transpose_conj()
    )

    pl2 = (1 - trace(outpostsel2 * ideal15to1)).real

    # Compute l2time, the speed at which level-3 rotations can be performed (t_{L2} in the 'paper')
    l2time = max(6 * dm2 / (nl2 / 2) / (1 - pfail2), dm3)

    # Define l2move, the effective width-dm3 region a level-2 state needs to traverse
    # before reaching the level-3 block, picking up additional storage errors
    lmove2 = 10 * dm3 + nl2 / 4 * (dx2 + 4 * dz2)

    # Step 1 of (15-to-1)x(15-to-1)x(8-to-CCZ) protocol applying rotations 1-2
    # Last operation: apply additional storage errors due to multi-patch measurements
    out3 = apply_rot(
        init4qubit,
        [z, one, one, z],
        pl2 + 0.5 * lmove2 * pm3,
        0.5 * lmove2 * pm3 + 0.5 * (3 * dx3 + dz3 + dm3) * dx3 / dm3 * pm3,
        mp.mpc(0),
    )
    out3 = apply_rot(
        out3,
        [one, one, one, -1 * z],
        pl2 + 0.5 * lmove2 * pm3,
        0.5 * lmove2 * pm3 + 0.5 * (dz3 + dm3) * dx3 / dm3 * pm3,
        mp.mpc(0),
    )
    out3 = storage_z_4(out3, 0.5 * (3 * dx3 + dz3 + dm3) * dm3 / dx3 * px3, 0, 0, 0)

    # Apply storage errors for l2time code cycles
    out3 = storage_x_4(out3, 0.5 * px3 * l2time, 0, 0, 0.5 * (dz3 / dx3) * px3 * l2time)
    out3 = storage_z_4(out3, 0.5 * px3 * l2time, 0, 0, 0.5 * (dx3 / dz3) * pz3 * l2time)

    # Step 2: apply rotations 3-4
    out3 = apply_rot(
        out3,
        [-1 * z, z, one, z],
        pl2 + 0.5 * lmove2 * pm3,
        0.5 * lmove2 * pm3 + 0.5 * (3 * dx3 + dz3 + dm3) * dx3 / dm3 * pm3,
        mp.mpc(0),
    )
    out3 = apply_rot(
        out3,
        [-1 * z, one, z, z],
        pl2 + 0.5 * lmove2 * pm3,
        0.5 * lmove2 * pm3 + 0.5 * (3 * dx3 + dz3 + dm3) * dx3 / dm3 * pm3,
        mp.mpc(0),
    )
    out3 = storage_z_4(
        out3,
        0.5 * ((3 * dx3 + dz3 + dm3) + (3 * dx3 + dz3 + dm3)) * dm3 / dx3 * px3,
        0.5 * (3 * dx3 + dz3 + dm3) * dm3 / dx3 * px3,
        0.5 * (3 * dx3 + dz3 + dm3) * dm3 / dx3 * px3,
        0,
    )

    # Apply storage errors for l2time code cycles
    out3 = storage_x_4(
        out3,
        0.5 * px3 * l2time,
        0.5 * px3 * l2time,
        0.5 * px3 * l2time,
        0.5 * (dz3 / dx3) * px3 * l2time,
    )
    out3 = storage_z_4(
        out3,
        0.5 * px3 * l2time,
        0.5 * px3 * l2time,
        0.5 * px3 * l2time,
        0.5 * (dx3 / dz3) * pz3 * l2time,
    )

    # Step 3: apply rotations 5-6
    out3 = apply_rot(
        out3,
        [z, z, z, z],
        pl2 + 0.5 * lmove2 * pm3,
        0.5 * lmove2 * pm3 + 0.5 * (3 * dx3 + dz3 + dm3) * dx3 / dm3 * pm3,
        mp.mpc(0),
    )
    out3 = apply_rot(
        out3,
        [one, -1 * z, z, z],
        pl2 + 0.5 * lmove2 * pm3,
        0.5 * lmove2 * pm3 + 0.5 * (2 * dx3 + dz3 + dm3) * dx3 / dm3 * pm3,
        mp.mpc(0),
    )
    out3 = storage_z_4(
        out3,
        0.5 * (3 * dx3 + dz3 + dm3) * dm3 / dx3 * px3,
        0.5 * ((3 * dx3 + dz3 + dm3) + (2 * dx3 + dz3 + dm3)) * dm3 / dx3 * px3,
        0.5 * ((3 * dx3 + dz3 + dm3) + (2 * dx3 + dz3 + dm3)) * dm3 / dx3 * px3,
        0,
    )

    # Apply storage errors for l2time code cycles
    # Qubit 1 is consumed as an output state: additional storage errors for dx3 code cycles
    out3 = storage_x_4(
        out3,
        0.5 * px3 * (dm3 + 2 * dx3),
        0.5 * px3 * l2time,
        0.5 * px3 * l2time,
        0.5 * (dz3 / dx3) * px3 * l2time,
    )
    out3 = storage_z_4(
        out3,
        0.5 * px3 * (dm3 + 2 * dx3),
        0.5 * px3 * l2time,
        0.5 * px3 * l2time,
        0.5 * (dx3 / dz3) * pz3 * l2time,
    )

    # Step 4: apply rotations 7-8
    out3 = apply_rot(
        out3,
        [one, z, one, z],
        pl2 + 0.5 * lmove2 * pm3,
        0.5 * lmove2 * pm3 + 0.5 * (3 * dx3 + dz3 + dm3) * dx3 / dm3 * pm3,
        mp.mpc(0),
    )
    out3 = apply_rot(
        out3,
        [one, one, z, z],
        pl2 + 0.5 * lmove2 * pm3,
        0.5 * lmove2 * pm3 + 0.5 * (dx3 + dz3 + dm3) * dx3 / dm3 * pm3,
        mp.mpc(0),
    )
    out3 = storage_z_4(
        out3,
        0,
        0.5 * (3 * dx3 + dz3 + dm3) * dm3 / dx3 * px3,
        0.5 * (dx3 + dz3 + dm3) * dm3 / dx3 * px3,
        0,
    )

    # Apply storage errors for l2time code cycles
    # Qubit 2 is consumed as an output state: additional storage errors for dx3 code cycles
    # Qubit 3 is consumed as an output state in the following step: additional storage errors for dx3 code cycles
    out3 = storage_x_4(
        out3,
        0,
        0.5 * px3 * (dm3 + 2 * dx3),
        0.5 * px3 * l2time + 0.5 * px3 * (dm3 + 2 * dx3),
        0.5 * (dz3 / dx3) * px3 * l2time,
    )
    out3 = storage_z_4(
        out3,
        0,
        0.5 * px3 * (dm3 + 2 * dx3),
        0.5 * px3 * l2time + 0.5 * px3 * (dm3 + 2 * dx3),
        0.5 * (dx3 / dz3) * pz3 * l2time,
    )

    return out3, l2time



def cost_of_three_level_8toccz(
    pphys: float | mpmath.mpf,
    dx: int,
    dz: int,
    dm: int,
    dx2: int,
    dz2: int,
    dm2: int,
    dx3: int,
    dz3: int,
    dm3: int,
    nl1: int,
    nl2: int,
) -> MagicStateFactory:
    """
    Generates the output-state density matrix of the (15-to-1)x(15-to-1)x(8-to-CCZ) protocol Calculates the output error and cost of the (15-to-1)x(15-to-1)x(8-to-CCZ) protocol with a physical error rate pphys, level-1 distances `dx`, `dz` and `dm`, level-2 distances `dx2, `dz2` and `dm2`, using `nl1` level-1 factories
    """

    pphys = mp.mpf(pphys)
    
    # Generate output state of (15-to-1)x(15-to-1)x(8-to-CCZ) protocol
    out3, l2time = three_level_8toccz_state(pphys, dx, dz, dm, dx2, dz2, dm2, dx3, dz3, dm3, nl1, nl2)

    # Compute level-3 failure probability as the probability to measure qubit 4 in the |+> state
    pfail3 = (1 - trace(kron(one, one, one, projx) * out3)).real

    # Compute the density matrix of the post-selected output state, i.e., after projecting qubit 4 into |+>
    outpostsel3 = (
        (1 / (1 - pfail3))
        * kron(one, one, one, projx)
        * out3
        * kron(one, one, one, projx).transpose_conj()
    )

    # Compute level-3 output error from the infidelity between the post-selected state and the ideal output state
    pout3 = (1 - trace(outpostsel3 * ideal8toCCZ)).real

    # Full-distance computation: determine full distance required for a 100-qubit / 10000-qubit computation
    def logerr1(d):
        return float(231 / (pout3 / 4) * d * plog(pphys, d) - 0.01)

    def logerr2(d):
        return float(20284 / (pout3 / 4) * d * plog(pphys, d) - 0.01)

    reqdist1 = int(2 * round(optimize.root(logerr1, 3, method="hybr").x[0] / 2) + 1)
    reqdist2 = int(2 * round(optimize.root(logerr2, 3, method="hybr").x[0] / 2) + 1)

    # Print output error, failure probability, space cost, time cost and space-time cost
    nqubits = 2 * int(
        (3 * dx3 + dz3) * 3 * dx3
        + nl2 * (dx2 + 4 * dz2) * dm3 / 2
        + 20 * dm3 * dm3
        + 2 * dx3 * dm3
    ) + 2 * nl2 * int(
        (3 * dx2 + dz2) * 3 * dx2
        + nl1 * ((dx + 4 * dz) * (3 * dx + dm2 / 2) + 2 * dm)
        + 20 * dm2 * dm2
        + 2 * dx2 * dm2
    )
    ncycles = 4 * l2time / (1 - pfail3)

    return MagicStateFactory(
        name=f"(15-to-1)x(15-to-1)x(8-to-CCZ) with pphys={float(pphys)}, dx={dx}, dz={dz}, dm={dm}, dx2={dx2}, dz2={dz2}, dm2={dm2}, dx3={dx3}, dz3={dz3}, dm3={dm3}, nl1={nl1}, nl2={nl2}",
        distilled_magic_state_error_rate=float(pout3),
        qubits=nqubits,
        distillation_time_in_cycles=float(ncycles),
        n_t_gates_produced_per_distillation=1,
    )
