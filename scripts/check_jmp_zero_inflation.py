#!/usr/bin/env python3
"""Independent SciPy/mpmath oracle through exported C ABI and loaded SQL.

Uses direct library distributions and high precision gamma/beta likelihoods,
not production recurrence/automatic differentiation. Native HC0 is explicit;
matching JMP's undocumented finite-sample sandwich details is not claimed.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import mpmath as mp
import numpy as np
import scipy
from scipy import optimize, special, stats

from check_jmp_censored import Abi
from check_jmp_continuation import Runner, expression

ROOT = Path(__file__).resolve().parents[1]
NAMES = ["zi_binomial", "zi_beta_binomial", "zi_poisson", "zi_negative_binomial", "zi_gamma"]
DISP = {1, 3, 4}


def rv(family, eta, dispersion, trials):
    mu = special.expit(eta) if family < 2 else np.exp(eta)
    if family == 0:
        return stats.binom(trials, mu)
    if family == 1:
        return stats.betabinom(trials, mu*(1/dispersion-1), (1-mu)*(1/dispersion-1))
    if family == 2:
        return stats.poisson(mu)
    if family == 3:
        return stats.nbinom(1/dispersion, 1/(1+mu*dispersion))
    return stats.gamma(mu/dispersion, scale=dispersion)


def fixtures():
    for family in range(5):
        rows = []
        for x in [-1., 0., 1.]:
            eta = (-0.4 if family < 2 else 0.5)+0.3*x
            dispersion = 0.2 if family == 1 else 0.5
            trials = 6 if family < 2 else 1
            if family == 4:
                rows.append([x, 0., 1., 7.])
                rows.extend([[x, math.exp(eta)*v, 1., w] for v, w in zip([.4, .7, 1., 1.6, 2.1], [2, 3, 5, 3, 2])])
            else:
                r = rv(family, eta, dispersion, trials)
                for y in range(7 if family < 2 else 30):
                    frequency = round(10000*((.3 if y == 0 else 0)+.7*r.pmf(y)))
                    if frequency:
                        rows.append([x, float(y), float(trials), float(frequency)])
        yield family, rows


def physical(f, t):
    beta = np.asarray(t[:-2] if f in DISP else t[:-1])
    d = (special.expit(t[-2]) if f == 1 else math.exp(t[-2])) if f in DISP else 0.
    return beta, d, special.expit(t[-1])


def scipy_reference(f, rows):
    a = np.asarray(rows)
    x = np.column_stack([np.ones(len(a)), a[:, :-3]])
    y, trials, w = a[:, -3:].T
    def fun(t):
        if np.any(np.abs(t) > 25):
            return 1e50
        beta, d, pi = physical(f, t)
        r = rv(f, x@beta, d, trials)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            l = r.logpdf(y) if f == 4 else r.logpmf(y)
            l = math.log1p(-pi)+l
            l[y == 0] = math.log(pi) if f == 4 else np.logaddexp(math.log(pi), l[y == 0])
        value = -float(w@l)/w.sum()
        return value if math.isfinite(value) else 1e50
    candidates = []
    for inflation in [.2, .6]:
        for dispersion in ([.15, .65] if f in DISP else [0.]):
            start = [0., 0.]+([special.logit(dispersion) if f == 1 else math.log(dispersion)] if f in DISP else [])+[special.logit(inflation)]
            first = optimize.minimize(fun, start, method="BFGS", options={"gtol": 1e-9, "maxiter": 1000})
            candidates.append(first)
    best = min(candidates, key=lambda r: r.fun)
    beta, d, pi = physical(f, best.x)
    return np.r_[beta, [d] if f in DISP else [], pi], best


def mp_reference(f, rows, start):
    mp.mp.dps = 60
    a = [[mp.mpf(v) for v in row] for row in rows]
    p = len(rows[0])-3
    def rowlog(t, row):
        eta = t[0]+sum(t[j+1]*row[j] for j in range(p))
        mu = 1/(1+mp.exp(-eta)) if f < 2 else mp.exp(eta)
        d = t[-2] if f in DISP else mp.mpf(0)
        pi = t[-1]
        y, n = row[-3:-1]
        if f == 4:
            if y == 0:
                return mp.log(pi)
            density = (y/d)**(mu/d-1)*mp.exp(-y/d)/(d*mp.gamma(mu/d))
            return mp.log1p(-pi)+mp.log(density)
        choose = mp.binomial(n, y) if f < 2 else 1
        if f == 0:
            density = choose*mu**y*(1-mu)**(n-y)
        elif f == 1:
            aa, bb = mu*(1/d-1), (1-mu)*(1/d-1)
            density = choose*mp.beta(y+aa, n-y+bb)/mp.beta(aa, bb)
        elif f == 2:
            density = mp.exp(-mu)*mu**y/mp.gamma(y+1)
        else:
            density = mp.gamma(y+1/d)/(mp.gamma(1/d)*mp.gamma(y+1))*(1/(1+mu*d))**(1/d)*(mu*d/(1+mu*d))**y
        return mp.log((pi if y == 0 else 0)+(1-pi)*density)
    def likelihood(*theta):
        return sum(row[-1]*rowlog(theta, row) for row in a)
    size = len(start)
    def gradient(*theta):
        return tuple(mp.diff(lambda v: likelihood(*(v if j == i else t for j, t in enumerate(theta))), theta[i]) for i in range(size))
    theta = mp.findroot(gradient, tuple(mp.mpf(v) for v in start), tol=mp.mpf("1e-42"), maxsteps=40)
    orders = lambda i, j: tuple(int(k == i)+int(k == j) for k in range(size))
    h = mp.matrix([[-mp.diff(likelihood, tuple(theta), orders(i, j)) for j in range(size)] for i in range(size)])
    model = h**-1
    meat = mp.matrix(size)
    for row in a:
        g = mp.matrix([mp.diff(lambda v: rowlog(tuple(v if j == i else t for j, t in enumerate(theta)), row), theta[i]) for i in range(size)])
        meat += row[-1]*(g*g.T)
    hc0 = model*meat*model
    tolist = lambda m: [[float(m[i, j]) for j in range(m.cols)] for i in range(m.rows)]
    result = {"theta": [float(v) for v in theta], "log_likelihood": float(likelihood(*theta)), "model_based": tolist(model), "sandwich_HC0": tolist(hc0), "score_max": float(max(abs(g) for g in gradient(*theta))), "mpmath_dps": mp.mp.dps}
    assert np.linalg.eigvalsh(np.array(result["model_based"])).min() > 0
    assert result["score_max"] < 1e-25
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--duckdb")
    parser.add_argument("--extension", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    abi = Abi(args.library.resolve())
    boundaries = {"exported_C_ABI": abi}
    if args.duckdb:
        assert args.extension
        sql = Runner(args.duckdb, args.extension)
        boundaries["loaded_DuckDB_SQL"] = lambda op, rows, parameters: sql.one(expression(op, rows, parameters))
    checks = 0
    errors = []
    def same(got, want, name, tol=3e-5, relative=False):
        nonlocal checks
        g, w = np.asarray(got, dtype=float), np.asarray(want, dtype=float)
        assert g.shape == w.shape, (name, g.shape, w.shape)
        assert np.isfinite(g).all(), (name, g)
        denominator = np.maximum(np.abs(w), 1e-300) if relative else 1+np.abs(w)
        scaled = float(np.max(np.abs(g-w)/denominator))
        errors.append(scaled)
        assert scaled <= tol, (name, scaled, got, want)
        checks += g.size
    cases = []
    for f, rows in fixtures():
        start, optimizer = scipy_reference(f, rows)
        reference = mp_reference(f, rows, start)
        theta = np.asarray(reference["theta"])
        beta, disp, pi = theta[:2], (theta[-2] if f in DISP else 0.), theta[-1]
        case = {"family": NAMES[f], "data": rows, "reference": reference, "scipy_initial_success": bool(optimizer.success), "boundaries": {}}
        for layer, run in boundaries.items():
            fits = []
            for kind, key in enumerate(["model_based", "sandwich_HC0"]):
                got = run(NAMES[f], rows, [600., 1e-9, .95, float(kind)])
                assert got["converged"] and got["inference_valid"] and not got["at_boundary"], (layer, NAMES[f], got)
                same(got["coefficients"], beta, layer+".coefficients", 3e-6)
                same(got["pi"], pi, layer+".pi", 3e-6)
                same(got["dispersion"], disp, layer+".dispersion", 3e-6)
                same(got["log_likelihood"], reference["log_likelihood"], layer+".likelihood", 1e-10)
                same(got["covariance"], reference[key], layer+".covariance", 3e-6)
                queries = [[x, 6. if f < 2 else 1., 1., .8] for x in [-.3, .6]]
                pred = run("zi_predict", queries, got["predict_parameters"])
                cov = np.asarray(reference[key])
                for i, (x, n, y, u) in enumerate(queries):
                    r = rv(f, beta[0]+beta[1]*x, disp, n)
                    base_mean, base_var = r.mean(), r.var()
                    mean = (1-pi)*base_mean
                    zero = pi if f == 4 else pi+(1-pi)*r.pmf(0)
                    same(pred["mixture_mean"][i], mean, layer+".mean", 3e-6)
                    same(pred["mixture_variance"][i], (1-pi)*base_var+pi*(1-pi)*base_mean**2, layer+".variance", 3e-6)
                    same(pred["zero_probability"][i], zero, layer+".zero", 3e-6)
                    same(pred["cdf"][i], pi+(1-pi)*r.cdf(y), layer+".cdf", 3e-6)
                    same(pred["survival"][i], (1-pi)*r.sf(y), layer+".sf", 3e-6)
                    same(pred["quantile"][i], 0. if u <= pi else r.ppf((u-pi)/(1-pi)), layer+".quantile", 3e-6)
                    field = "query_positive_density" if f == 4 else "query_probability_mass"
                    same(pred[field][i], (1-pi)*(r.pdf(y) if f == 4 else r.pmf(y)), layer+".query", 3e-6)
                    derivative = base_mean*(1-special.expit(beta[0]+beta[1]*x)) if f < 2 else base_mean
                    gradient = np.r_[[(1-pi)*derivative, (1-pi)*derivative*x], [0.] if f in DISP else [], -base_mean]
                    se = math.sqrt(float(gradient@cov@gradient))
                    same(pred["mean_standard_error"][i], se, layer+".mean_SE", 3e-6)
                    same(pred["mean_confidence_interval"][i], mean+stats.norm.ppf([.025,.975])*se, layer+".mean_CI", 3e-6)
                fits.append(got)
            case["boundaries"][layer] = fits
        cases.append(case)
        print("PASS", NAMES[f], "boundaries", len(boundaries), flush=True)
    def carrier(f, eta, dispersion, pi):
        size = 3 if f in DISP else 2
        boundary = pi == 0 or (f in {1,3} and dispersion == 0)
        return [1., float(f), 0., .95, 0., 0., float(boundary), eta, dispersion, pi]+[0.]*(size*size)
    edges = []
    for f, cdf, variance in [(0, 8.303379160420424e-87, 2.5490125531749533e-17), (1, 1.9555916413246998e-19, 5.0980251063499066e-17)]:
        edges.append(("saturated_logit_"+str(f), [[6.,1.,.5]], carrier(f,40.,.2 if f else 0.,0.), {"cdf":[cdf],"mixture_variance":[variance]}))
    edges += [
        ("tiny_gamma_quantile", [[1.,1.,1e-250]], carrier(4,0.,1.,0.), {"quantile":[1e-250]}),
        ("underflow_scaled_gamma", [[1.,float.fromhex('0x0.0000000000001p-1022'),.5]], carrier(4,math.log(1e153),1e154,0.), {"cdf":[1.9545095271844815e-48]}),
        ("right_poisson_tail", [[1.,100.,.5]], carrier(2,math.log(2),0.,.3), {"survival":[.7*stats.poisson(2).sf(100)]}),
    ]
    pi=.3;u=np.nextafter(pi,1.)
    edges.append(("just_above_atom",[[1.,1.,float(u)]],carrier(4,0.,1.,pi),{"quantile":[-math.log1p(-(u-pi)/(1-pi))]}))
    # mpmath75 independent integrals/loggamma; the large-shape CDF integral
    # changes variable x=a+sqrt(a)*t and integrates t=-30..query.
    for name,eta,y,cdf,sf,pdf in [
        ("tiny_shape_normalizer",-17.72753356339242,.1,.9999999635415212,3.6458478795160347e-8,1.809674773624767e-7),
        ("small_shape_normalizer",-13.815510557964274,2.078255658098977e-44,.9999,9.999999999998897e-5,4.811246374349484e37),
        ("large_shape_normalizer",16.012735135300492,9000000.,.5000443269220536,.49995567307794647,.00013298075890250745),
    ]:
        edges.append((name,[[1.,y,.999999999999]],carrier(4,eta,1.,0.),{"cdf":[cdf],"survival":[sf],"query_positive_density":[pdf]}))
    gamma_grid=[]
    for shape in [1e-6,.01,.1,1.,10.,1000.,100000.,9000000.]:
        eta=math.log(shape);actual_shape=math.exp(eta);dist=stats.gamma(a=actual_shape)
        for probability in [.0001,.5,.9999]:
            y=float(dist.ppf(probability))
            if y:
                gamma_grid.append((eta,y,float(dist.cdf(y)),float(dist.sf(y))))
    edge_results = {}
    invalid = []
    for f, op in enumerate(NAMES):
        n=6. if f<2 else 1.
        invalid.extend([
            (op,[[0.,n,10.]],[]),
            (op,[[1.,n,-1.],[2.,n,2.]],[]),
            (op,[[1.,n,.5],[2.,n,2.]],[]),
            (op,[[1.,0.,n,10.],[1.,1.,n,10.]],[]),
        ])
    invalid += [
        ("zi_predict",[[1.,.01,1e-5]],carrier(4,math.log(.01),1.,0.)),
        ("zi_predict",[[1.,1.,0.]],carrier(2,math.log(2),0.,.3)),
        ("zi_predict",[[2.,1.,.5]],carrier(2,math.log(2),0.,.3)),
        ("zi_predict",[[1.,1.,.5]],carrier(2,math.log(2),0.,.3)+[0.]),
    ]
    invalid_checks=0
    for layer, run in boundaries.items():
        edge_results[layer]={}
        for name, rows, parameters, expected in edges:
            got=run("zi_predict",rows,parameters)
            for key,want in expected.items():
                same(got[key],want,layer+"."+name+"."+key,tol=2e-9,relative=True)
            edge_results[layer][name]={key:got[key] for key in expected}
        for eta,y,cdf,sf in gamma_grid:
            got=run("zi_predict",[[1.,y,.999999999999]],carrier(4,eta,1.,0.))
            same(got["cdf"],[cdf],layer+".Gamma_grid_CDF",tol=2e-9,relative=True)
            same(got["survival"],[sf],layer+".Gamma_grid_SF",tol=2e-9,relative=True)
        for op,rows,parameters in invalid:
            if layer=="exported_C_ABI":
                assert abi.nullable(op,rows,parameters) is None,(op,rows)
            else:
                assert sql.query("SELECT ("+expression(op,rows,parameters)+") IS NULL AS valid;")[0]["valid"],(op,rows)
            invalid_checks+=1
        # Integer frequencies must be equivalent to actual row replication,
        # including unadjusted HC0 meat; squaring frequencies breaks this.
        grouped=[[0.,1.,10.],[1.,1.,3.],[2.,1.,4.],[3.,1.,2.],[4.,1.,1.]]
        expanded=[row[:-1]+[1.] for row in grouped for _ in range(int(row[-1]))]
        first=run("zi_poisson",grouped,[600.,1e-9,.95,1.])
        second=run("zi_poisson",expanded,[600.,1e-9,.95,1.])
        for key in ["coefficients","pi","covariance","log_likelihood"]:
            same(first[key],second[key],layer+".frequency_replication."+key,tol=2e-7)
    sources=["crates/analytics/src/lib.rs","crates/analytics/src/zero_inflation.rs","crates/analytics/tests/zero_inflation_contracts.rs","src/include/analytics_registry.inc"]
    result = {"schema_version": 1, "status": "PASS", "scope": "Five bounded natural-support shared-inflation ML families, model-based and HC0 uncertainty, portable prediction; not full JMP equivalence", "python": sys.version, "platform": platform.platform(), "scipy": scipy.__version__, "mpmath": mp.__version__, "boundaries": list(boundaries), "checks": checks, "invalid_NULL_checks":invalid_checks, "gamma_shape_grid_cases_per_boundary":len(gamma_grid), "maximum_scaled_error": max(errors), "library_sha256": hashlib.sha256(args.library.read_bytes()).hexdigest(), "source_sha256":{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}, "numerical_primary_references":["https://dlmf.nist.gov/5.7","https://dlmf.nist.gov/5.11","https://dlmf.nist.gov/8.7"], "edge_cases":edge_results, "cases": cases}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k != "cases"}, indent=2))


if __name__ == "__main__":
    main()
