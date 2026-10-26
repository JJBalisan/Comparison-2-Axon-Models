"""Inhibition and the coincidence window (issue #77): width and centre vs dt_inh.

    uv run scripts/inhibition_window.py --out-dir figures/inhibition

An IPSG (rise 0.28 ms, decay 1.5 ms, reversal -80 mV, 40 nS by default) is
timed `dt_inh` ms after the first EPSG of a pair; the second EPSG is `delay`
after it, so the window is measured on a signed delay (coincidence.window_edges)
with Myoga et al 2014's EPSG (decay 0.3 ms) at 3% above the lowest threshold
over delay. For each dt_inh the table gives the window's width and centre (the
optimum Myoga et al found shifting), for

- the 45-compartment model with KLT present (the default membrane);
- the same with the rest balanced (membrane(rebalance_rest=True)), the
  control for the next column;
- the same with KLT replaced by an equal leak at rest (the KLT conductance at
  the resting potential, with its gating at steady state, added to the leak),
  which removes the Kv1 history effect Roberts et al 2013 invoke without
  changing the resting conductance;
- the two-compartment model of TwoCpt.m, KLT present (it has no knob for
  taking KLT out).

The second block holds the excitation at 3% above the threshold *without*
inhibition, as a fixed-size input would be: where the IPSG stops that input
spiking altogether there is no window.

Targets: Myoga et al 2014 (slice, no KLT history effect in their model): AP
half-width 208 us with inhibition at -0.2 ms against 221 us without; optimum
shifted -73 us with inhibition leading by 0.1 ms and +41 us lagging by 0.6 ms.
Roberts et al 2013, Franken et al 2015 and Casarez et al 2026: narrowing without
a shift of the centre. In vivo IPSG kinetics are not in the literature: the
kinetics here are slice values.
"""  # noqa: E501 (printed verbatim as --help)

import argparse
import dataclasses
import json
import pathlib
from concurrent import futures
from typing import Any

import numpy as np

import msoaxon
from msoaxon import _fast, coincidence, constants, measure, multi, synaptic

EPSG = (0.1, 0.3)  # Myoga et al 2014
MARGIN = 0.03
DT_INH = (-0.6, -0.4, -0.2, -0.1, 0.0, 0.1, 0.2, 0.4, 0.6)
SIGNED = "the second EPSG's delay after the first"


def klt_as_leak(v0: float = -68.0) -> multi.Morphology:
    """Returns LUMPED with each compartment's KLT moved into its leak.

    The leak density gains g_klt * w^4 * z at v0 (w and z at their steady
    state there), so the membrane's conductance at rest is unchanged and
    the voltage-dependent KLT is gone.

    Args:
        v0: The resting potential the conductance is evaluated at [mV].

    Returns:
        A Morphology like multi.LUMPED with no KLT, and its own key.
    """
    base = multi.LUMPED
    v = np.full(base.n, float(v0))
    at_rest = base.g_klt * constants.winf(v) ** 4 * constants.zinf(v)
    return dataclasses.replace(
        base,
        key="lumped, KLT replaced by an equal leak",
        g_klt=np.zeros(base.n),
        g_lk=base.g_lk + at_rest,
    )


def configs(v0: float) -> dict[str, tuple[str, float, dict[str, Any]]]:
    """Returns the model configurations: name -> (model, v0, model_kw)."""
    leak = klt_as_leak(v0)
    return {
        "45-compartment, KLT present": ("multi", v0, {}),
        "45-compartment, KLT present, rest balanced": (
            "multi",
            v0,
            {"mem": multi.membrane(v0, rebalance_rest=True)},
        ),
        "45-compartment, KLT replaced by equal leak": (
            "multi",
            v0,
            {
                "morph": leak,
                "mem": multi.membrane(v0, rebalance_rest=True, morph=leak),
            },
        ),
        "2-cpt, TwoCpt.m, KLT present": ("two", v0, {}),
    }


def edges(
    model: str,
    v0: float,
    model_kw: dict[str, Any],
    dt: float | None,
    g_inh: float,
    tol: float,
    delay_tol: float,
    reference: tuple[float, float] | None = None,
) -> tuple[float, float, float, float]:
    """Runs window_edges for one configuration and dt_inh (a pool job).

    Returns:
        The WindowEdges as a tuple (lo, hi, best, th0), delays in ms.
    """
    inh = None if dt is None else synaptic.Inhibition(dt, g_inh)
    src = synaptic.PairSource(tau=EPSG, inh=inh)
    return tuple(
        coincidence.window_edges(
            model,
            MARGIN,
            v0=v0,
            epsg_tau=EPSG,
            model_kw=model_kw,
            rel_tol=tol,
            delay_tol=delay_tol,
            reference=reference,
            events=src,
        )
    )


def noisy_check(
    a: argparse.Namespace,
    pool: futures.Executor,
    cfgs: dict[str, tuple[str, float, dict[str, Any]]],
    jobs: dict[tuple[str, float | None], futures.Future],
    dts: tuple[float, ...],
) -> dict[str, float]:
    """Checks window_edges against noisy trials (two-compartment, one dt_inh).

    Trials jitter each EPSG's onset by 5 us and the amplitude by 1% (as in
    scripts/coincidence_window.py), at 3% above the lowest threshold; the
    width and centre of the resulting probability curve are compared with
    window_edges'.
    """
    name = "2-cpt, TwoCpt.m, KLT present"
    dt = -0.2 if -0.2 in dts else dts[0]
    model, v0, kw = cfgs[name]
    lo, hi, best, th0 = jobs[name, dt].result()
    delays = np.round(np.arange(-0.4, 0.4001, 0.04 if not a.quick else 0.1), 3)
    prob = coincidence.probability_trials(
        model,
        delays,
        (1 + MARGIN) * th0,
        n_trials=a.trials if not a.quick else 4,
        amp_cv=0.01,
        jitter=0.005,
        v0=v0,
        epsg_tau=EPSG,
        model_kw=kw,
        executor=pool,
        dt_inh=dt,
        inh=synaptic.Inhibition(dt, a.g_inh),
    )
    s = measure.curve_summary(delays, prob)
    print(
        f"\nnoisy-trial check ({name}, dt_inh {dt:+g} ms): peak probability "
        f"{s['peak']:.2f}, half-width {s['width'] * 1e3:.0f} us, centre "
        f"{s['centre'] * 1e3:+.0f} us; window_edges: "
        f"{(hi - lo) * 1e3:.0f} us, {(hi + lo) / 2 * 1e3:+.0f} us"
    )
    return {"dt_inh": dt, **s, "window_width": hi - lo, "window_lo": lo}


def main() -> None:
    """Parses the command line and runs the inhibition sweep."""
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--out-dir", default="figures/inhibition")
    ap.add_argument("--g-inh", type=float, default=synaptic.IPSG_G)
    ap.add_argument(
        "--mode",
        choices=("fast", "native", "exact"),
        help="mso_axon's mode for the 45-compartment model",
    )
    ap.add_argument("--trials", type=int, default=30)
    ap.add_argument(
        "--quick",
        action="store_true",
        help="few dt_inh values and 1%% tolerances: a smoke test, not results",
    )
    a = ap.parse_args()
    with msoaxon.process_pool() as pool:
        run(a, pool)


def run(a: argparse.Namespace, pool: futures.Executor) -> None:
    """Runs the sweep, prints the tables and saves the JSON.

    Args:
        a: The parsed command-line arguments.
        pool: The process pool every parallel call shares.
    """
    out = pathlib.Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tol = 1e-2 if a.quick else 1e-4
    delay_tol = 2e-3 if a.quick else 1e-3  # [ms]
    dts = (-0.2, 0.2) if a.quick else DT_INH
    mode = _fast.effective(a.mode)
    cfgs = configs(-68.0)
    if a.mode:  # the mode applies to the 45-compartment model only
        for model, _, kw in cfgs.values():
            if model == "multi":
                kw["mode"] = a.mode
    print(
        f"IPSG {a.g_inh:g} nS, rise/decay {synaptic.IPSG_TAU}, reversal "
        f"{synaptic.IPSG_E_REV:g} mV; EPSG {EPSG}, margin {MARGIN:.0%}; "
        f"45-compartment mode {mode!r}; delay = {SIGNED}"
    )
    jobs = {
        (name, dt): pool.submit(
            edges, model, v0, kw, dt, a.g_inh, tol, delay_tol
        )
        for name, (model, v0, kw) in cfgs.items()
        for dt in (None, *dts)
    }
    res: dict[str, Any] = {"mode": mode, "g_inh": a.g_inh, "margin": MARGIN}
    first = {name: jobs[name, None].result() for name in cfgs}
    # the same input held at 3% above the no-inhibition threshold
    fixed = {
        (name, dt): pool.submit(
            edges,
            model,
            v0,
            kw,
            dt,
            a.g_inh,
            tol,
            delay_tol,
            (jobs[name, dt].result()[2], first[name][3]),
        )
        for name, (model, v0, kw) in cfgs.items()
        for dt in dts
    }
    print(
        "\nExcitation 3% above the lowest threshold with the IPSG present"
        "\nmodel | dt_inh [ms] | lowest threshold [nS] | width [us] | "
        "centre [us] (shift vs no inhibition)"
    )
    for name in cfgs:
        lo, hi, best, th0 = first[name]
        print(
            f"{name} | none | {th0:.2f} | {(hi - lo) * 1e3:.0f} | "
            f"{(hi + lo) / 2 * 1e3:+.0f}",
            flush=True,
        )
        res[name] = {
            "none": dict(
                zip(("lo", "hi", "best", "th0"), first[name], strict=True)
            )
        }
        for dt in dts:
            lo, hi, best, th0 = jobs[name, dt].result()
            print(
                f"{name} | {dt:+g} | {th0:.2f} | {(hi - lo) * 1e3:.0f} | "
                f"{(hi + lo) / 2 * 1e3:+.0f}",
                flush=True,
            )
            res[name][f"{dt:+g}"] = dict(
                zip(
                    ("lo", "hi", "best", "th0"),
                    (lo, hi, best, th0),
                    strict=True,
                )
            )
    print(
        "\nExcitation held at 3% above the threshold without inhibition"
        "\nmodel | dt_inh [ms] | width [us] | centre [us]"
    )
    for name in cfgs:
        for dt in dts:
            lo, hi, _, _ = fixed[name, dt].result()
            print(
                f"{name} | {dt:+g} | "
                + (
                    f"{(hi - lo) * 1e3:.0f} | {(hi + lo) / 2 * 1e3:+.0f}"
                    if np.isfinite(hi - lo)
                    else "no window | -"
                ),
                flush=True,
            )
            res[name][f"{dt:+g}"]["fixed_excitation"] = {"lo": lo, "hi": hi}
    res["noisy_check"] = noisy_check(a, pool, cfgs, jobs, dts)
    (out / "inhibition_window.json").write_text(json.dumps(res, indent=1))
    print(f"\nsaved {out}/inhibition_window.json")


if __name__ == "__main__":
    main()
