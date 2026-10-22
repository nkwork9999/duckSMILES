#!/usr/bin/env python3
"""Independent JMP19.2 censored-family SciPy oracle and public ABI/native SQL checks.

The oracle uses SciPy distributions and optimize.minimize (not Rust functions),
with numerical Hessian independently differenced in original coefficient units.
No reference-only run is presented as production or SQL verification.
"""
import argparse
import ctypes
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
import scipy
from scipy import optimize, stats
from check_jmp_continuation import strict_json

FAMILIES = ["normal", "exponential", "gamma", "weibull", "lognormal", "beta"]
PRIMARY = [
    "https://www.jmp.com/support/help/en/19.2/jmp/launch-the-generalized-regression-personality.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/specify-a-distribution.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/statistical-details-for-distributions.shtml",
]


def distribution(family, eta, log_sigma):
    sigma = np.exp(log_sigma)
    if family == "normal":
        return stats.norm(loc=eta, scale=sigma)
    if family == "exponential":
        return stats.expon(scale=np.exp(eta))
    if family == "gamma":
        # Official JMP19.2 density: shape mu/sigma and scale sigma.
        return stats.gamma(a=np.exp(eta)/sigma, scale=sigma)
    if family == "weibull":
        return stats.weibull_min(c=1/sigma, scale=np.exp(eta))
    if family == "lognormal":
        return stats.lognorm(s=sigma, scale=np.exp(eta))
    mu = scipy.special.expit(eta)
    return stats.beta(a=mu/sigma, b=(1-mu)/sigma)


def log_difference(a, b):
    if np.isneginf(b):
        return a
    if a <= b:
        return -np.inf
    return a + np.log(-np.expm1(b-a))


def nll(family, data, theta):
    data = np.asarray(data)
    p = data.shape[1]-4
    eta = theta[0] + data[:, :p] @ theta[1:p+1]
    log_sigma = 0.0 if family == "exponential" else theta[-1]
    rv = distribution(family, eta, log_sigma)
    lo, hi, kinds, weights = data[:, p:].T
    result = np.zeros(len(data))
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for i, kind in enumerate(kinds.astype(int)):
            if not weights[i]:
                continue
            if kind == 0:
                result[i] = rv.logpdf(lo)[i]
            elif kind == 1:
                result[i] = rv.logcdf(hi)[i]
            elif kind == 2:
                result[i] = rv.logsf(lo)[i]
            else:
                a, b = rv.logcdf(hi)[i], rv.logcdf(lo)[i]
                c, d = rv.logsf(lo)[i], rv.logsf(hi)[i]
                result[i] = log_difference(a, b) if b < -math.log(2) else log_difference(c, d)
    v = -float(result @ weights)
    return v if math.isfinite(v) else 1e100


def numerical_hessian(fun, x, step=4e-4):
    x = np.array(x, dtype=float)
    h = np.empty((len(x), len(x)))
    f0 = fun(x)
    for i in range(len(x)):
        for j in range(i+1):
            estimates = []
            for divisor in [1, 2]:
                a = step*(1+abs(x[i]))/divisor
                b = step*(1+abs(x[j]))/divisor
                ei, ej = np.eye(len(x))[i]*a, np.eye(len(x))[j]*b
                if i == j:
                    value = (fun(x+ei)-2*f0+fun(x-ei))/a**2
                else:
                    value = (fun(x+ei+ej)-fun(x+ei-ej)-fun(x-ei+ej)+fun(x-ei-ej))/(4*a*b)
                estimates.append(value)
            h[i, j] = h[j, i] = (4*estimates[1]-estimates[0])/3
    return h


def reference_fit(family, data):
    a = np.array(data, dtype=float)
    p = a.shape[1]-4
    proxy = np.where(a[:, p+2] == 1, a[:, p+1], a[:, p])
    proxy = np.where(a[:, p+2] == 3, (a[:, p]+a[:, p+1])/2, proxy)
    if family in ["weibull", "lognormal"]:
        proxy = np.log(np.maximum(proxy, 1e-8))
    intercept = float(np.average(proxy, weights=a[:, -1]))
    if family in ["exponential", "gamma"]:
        intercept = math.log(intercept)
    if family == "beta":
        intercept = math.log(intercept/(1-intercept))
    initial = [intercept] + [0.0]*p + ([] if family == "exponential" else [-0.5])
    fun = lambda theta: nll(family, a, theta)
    candidates = []
    for delta in [-0.7, 0.0, 0.7]:
        start = np.array(initial)
        if family != "exponential":
            start[-1] += delta
        # Powell search is deliberately different from production BFGS/FD.
        first = optimize.minimize(fun, start, method="Powell", options={"xtol": 1e-11, "ftol": 1e-13, "maxiter": 3000})
        second = optimize.minimize(fun, first.x, method="BFGS", jac="3-point", options={"gtol": 1e-7, "maxiter": 1000})
        candidates.extend([first, second])
    best = min(candidates, key=lambda r: r.fun)
    hessian = numerical_hessian(fun, best.x)
    eigenvalues = np.linalg.eigvalsh(hessian)
    assert np.min(eigenvalues)>0, (family, best, eigenvalues)
    cov = np.linalg.inv(hessian)
    return {
        "theta": best.x.tolist(), "log_likelihood": -float(best.fun),
        "standard_errors": np.sqrt(np.diag(cov)).tolist(), "covariance": cov.tolist(),
        "optimizer_success": bool(best.success), "optimizer_message": str(best.message),
        "observed_information_eigenvalues": eigenvalues.tolist(),
    }


def fixtures():
    # Mixed finite observations, all four censor kinds, two independent covariates.
    x = np.array([-1, -.8, -.5, -.2, 0, .2, .4, .6, .8, 1, 1.2, 1.4], dtype=float)
    z = np.array([.2, -.5, .8, -.8, .5, -.2, .9, -.9, .3, -.3, .6, -.6])
    noise = np.array([-.8, .4, -.2, .8, -.6, .3, -.3, .7, -.5, .2, -.1, .6])
    for family in FAMILIES:
        if family == "normal":
            y = .7+.5*x-.3*z+noise
        elif family == "beta":
            y = scipy.special.expit(-.3+.5*x-.4*z+noise)
        else:
            y = np.exp(.4+.35*x-.25*z+noise)
        rows = []
        for i, (a, b, v) in enumerate(zip(x, z, y)):
            if i < 7:
                lo, hi, kind = v, v, 0
            elif i < 9:
                lo, hi, kind = 0., v, 1
            elif i < 11:
                lo, hi, kind = v, 0., 2
            else:
                lo, hi, kind = (v-.05 if family == "beta" else v*.9), (v+.05 if family == "beta" else v*1.15), 3
            rows.append([float(a), float(b), float(lo), float(hi), kind, 2 if i%4==0 else 1])
        yield family, rows


class Abi:
    def __init__(self, path):
        self.lib = ctypes.CDLL(str(path))
        self.run = self.lib.ds_analytics_run
        self.run.argtypes = [ctypes.c_void_p,ctypes.c_size_t,ctypes.c_void_p,ctypes.c_size_t,ctypes.c_size_t,ctypes.c_void_p,ctypes.c_size_t,ctypes.POINTER(ctypes.c_size_t)]
        self.run.restype = ctypes.c_void_p
        self.free = self.lib.ds_analytics_free
        self.free.argtypes = [ctypes.c_void_p,ctypes.c_size_t]
    def nullable(self, operation, data, args):
        op = operation.encode()
        a = np.ascontiguousarray(data,dtype=np.float64)
        b = np.ascontiguousarray(args,dtype=np.float64)
        n = ctypes.c_size_t()
        ptr = self.run(op,len(op),a.ctypes.data,*a.shape,b.ctypes.data,len(b),ctypes.byref(n))
        if not ptr:
            assert n.value == 0, (operation, "NULL buffer with nonzero length", n.value)
            return None
        try:
            return strict_json(ctypes.string_at(ptr,n.value))
        finally:
            self.free(ptr,n.value)
    def __call__(self, operation, data, args):
        result = self.nullable(operation, data, args)
        assert result is not None, f"production ABI rejected {operation}: {data} {args}"
        return result


def confirmed_boundary_regressions(boundaries, nullable_boundaries):
    """Previously RED Rust cases, independently asserted at both compiled boundaries."""
    fixture_path = Path(__file__).resolve().parents[1]/"tasks/evidence/jmp-censored-gamma-tail-review.json"
    fixture = json.loads(fixture_path.read_text())
    gamma_data = [case["input_data"][0] for case in fixture["cases"]]
    gamma = {"name": "gamma_actual_small_upper_tail", "data": gamma_data,
             "predict_parameters": fixture["predict_parameters"],
             "independent_fixture": {"path": str(fixture_path.relative_to(fixture_path.parents[2])),
                                     "sha256": hashlib.sha256(fixture_path.read_bytes()).hexdigest(),
                                     "method": fixture["reference_method"],
                                     "environment": fixture["reference_environment"]},
             "results": {}}
    interval_cases = [
        ("exponential_overflowing_mean_ci", [[0., .1]], [1.,1.,0.,.95,0.,1.,1.,709.,0.,1.]),
        ("normal_overflowing_mean_ci", [[0., .5]], [1.,0.,0.,.95,1e308,1e308,1.,0.,0.,1.,0.,0.,1.]),
    ]
    intervals = [{"name": name, "data": data, "predict_parameters": params, "results": {}}
                 for name, data, params in interval_cases]
    asymmetric = {"name": "huge_asymmetric_covariance_rejected", "data": [[0., .5]],
                  "predict_parameters": [1.,0.,0.,.95,0.,1.,1.,0.,0.,1.7e308,1.1e308,1e308,1.7e308],
                  "results": {}}
    checks = 0
    for layer, call in boundaries.items():
        actual = call("censored_predict", gamma_data, fixture["predict_parameters"])
        for field, oracle, rtol, atol in [("cdf", "cdf", 0., 1e-15),
                                           ("survival", "survival", 1e-7, 0.),
                                           ("density", "density", 1e-7, 0.),
                                           ("quantiles", "quantile", 0., 1e-8)]:
            expected = [float(case["mpmath80dps"][oracle]) for case in fixture["cases"]]
            # Rust may emit a finite large float as an integer JSON token. Explicit
            # binary64 conversion avoids NumPy object arrays of arbitrary-size ints.
            np.testing.assert_allclose(np.asarray(actual[field], dtype=float), expected, rtol=rtol, atol=atol)
            checks += len(expected)
        gamma["results"][layer] = actual
        for case in intervals:
            actual = call("censored_predict", case["data"], case["predict_parameters"])
            assert len(actual["means"]) == 1, (layer, case["name"], actual)
            assert len(actual["mean_standard_errors"]) == 1, (layer, case["name"], actual)
            mean, error = actual["means"][0], actual["mean_standard_errors"][0]
            assert math.isfinite(mean) and mean > 0, (layer, case["name"], mean)
            assert math.isfinite(error) and error > 0, (layer, case["name"], error)
            assert actual["mean_confidence_intervals"] == [None], (layer, case["name"], actual)
            checks += 5
            case["results"][layer] = actual
        actual = nullable_boundaries[layer]("censored_predict", asymmetric["data"], asymmetric["predict_parameters"])
        assert actual is None, (layer, "accepted materially asymmetric covariance", actual)
        asymmetric["results"][layer] = {"returned_null": True}
        checks += 1
    return [gamma, *intervals, asymmetric], checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library",type=Path)
    parser.add_argument("--duckdb")
    parser.add_argument("--extension",type=Path)
    parser.add_argument("--reference-only",action="store_true")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if not (args.library or args.duckdb or args.reference_only):
        parser.error("select --library, --duckdb/--extension, or explicitly --reference-only")
    if bool(args.duckdb)!=bool(args.extension):
        parser.error("--duckdb and --extension must be supplied together")
    boundaries={}
    nullable_boundaries={}
    if args.library:
        abi=Abi(args.library)
        boundaries["public_c_abi"]=abi
        nullable_boundaries["public_c_abi"]=abi.nullable
    if args.duckdb:
        from check_jmp_continuation import Runner, expression
        runner=Runner(args.duckdb,args.extension)
        boundaries["native_duckdb_sql"]=lambda op,d,a:runner.one(expression(op,d,a))
        def sql_nullable(op, data, params):
            rows=runner.query("SELECT "+expression(op,data,params)+" AS result;")
            assert len(rows)==1, rows
            value=rows[0]["result"]
            if value is None:
                return None
            runner.decoded+=1
            return strict_json(value)
        nullable_boundaries["native_duckdb_sql"]=sql_nullable
    fit_args=[600,1e-7,.95]
    output={"schema":"jmp19.2-censored-independent-reference-v1","primary_sources":PRIMARY,
            "fit_arguments":fit_args,
            "fixture_precision":"Uses documented default optimizer gradient tolerance1e-7. Independent numerical reference acceptance tolerances are unchanged; returned fits additionally require an independent SciPy central-difference average-NLL score <=2e-7.",
            "reference_versions":{"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__},
            "platform":platform.platform(),"source_revision":subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
            "verification_boundaries":list(boundaries),"reference_only":not bool(boundaries),"cases":[],"assertions":0}
    for family,data in fixtures():
        reference=reference_fit(family,data)
        case={"family":family,"data":data,"independent_reference":reference,"results":{}}
        for layer,call in boundaries.items():
            result=call("censored_"+family,data,fit_args)
            assert result["converged"] and result["inference_valid"], (layer,family,result)
            actual_theta=np.array(result["coefficients"]+([] if family=="exponential" else [result["log_sigma"]]))
            total_frequency=sum(row[-1] for row in data)
            score=scipy.optimize._numdiff.approx_derivative(lambda t:nll(family,data,t)/total_frequency,actual_theta,method="3-point")
            score_max=float(np.max(np.abs(score)))
            assert score_max<=2e-7, (layer,family,"independent average-NLL score",score_max)
            output["assertions"]+=2
            np.testing.assert_allclose(result["coefficients"],reference["theta"][:3],rtol=2e-4,atol=3e-5)
            if family!="exponential":
                np.testing.assert_allclose(result["log_sigma"],reference["theta"][-1],rtol=2e-4,atol=3e-5)
            np.testing.assert_allclose(result["log_likelihood"],reference["log_likelihood"],rtol=1e-8,atol=1e-7)
            np.testing.assert_allclose(result["standard_errors"],reference["standard_errors"],rtol=5e-4,atol=2e-5)
            np.testing.assert_allclose(result["covariance"],reference["covariance"],rtol=1e-3,atol=3e-5)
            z95=stats.norm.ppf(.975)
            theta_ref=np.array(reference["theta"])
            se_ref=np.array(reference["standard_errors"])
            expected_ci=np.column_stack([theta_ref-z95*se_ref,theta_ref+z95*se_ref])
            np.testing.assert_allclose(result["confidence_intervals"],expected_ci,rtol=1e-3,atol=5e-5)
            output["assertions"]+=len(theta_ref)**2+2*len(theta_ref)
            query=[[0,0,.45 if family=="beta" else 1.,.1],[.5,-.3,.65 if family=="beta" else 2.,.5],[1,.2,.3 if family=="beta" else .5,.9]]
            pred=call("censored_predict",query,result["predict_parameters"])
            xx=np.array(query);theta=np.array(reference["theta"])
            eta=theta[0]+xx[:,:2]@theta[1:3]
            rv=distribution(family,eta,0 if family=="exponential" else theta[-1])
            for key,expected in [("means",rv.mean()),("cdf",rv.cdf(xx[:,2])),("survival",rv.sf(xx[:,2])),("density",rv.pdf(xx[:,2])),("quantiles",rv.ppf(xx[:,3]))]:
                np.testing.assert_allclose(pred[key],expected,rtol=3e-4,atol=3e-5)
                output["assertions"]+=len(expected)
            # Joint observed-information delta mean variance, calculated on raw coefficients.
            def fmean(t):
                return distribution(family,t[0]+xx[:,:2]@t[1:3],0 if family=="exponential" else t[-1]).mean()
            jac=scipy.optimize._numdiff.approx_derivative(fmean,theta,method="3-point")
            covariance=np.array(reference["covariance"])
            wanted=np.sqrt(np.einsum("ij,jk,ik->i",jac,covariance,jac))
            np.testing.assert_allclose(pred["mean_standard_errors"],wanted,rtol=1e-3,atol=3e-5)
            expected_mean_ci=np.column_stack([rv.mean()-z95*wanted,rv.mean()+z95*wanted])
            np.testing.assert_allclose(pred["mean_confidence_intervals"],expected_mean_ci,rtol=2e-3,atol=5e-5)
            def fquantile(t):
                return distribution(family,t[0]+xx[:,:2]@t[1:3],0 if family=="exponential" else t[-1]).ppf(xx[:,3])
            qjac=scipy.optimize._numdiff.approx_derivative(fquantile,theta,method="3-point")
            wanted_qse=np.sqrt(np.einsum("ij,jk,ik->i",qjac,covariance,qjac))
            np.testing.assert_allclose(pred["quantile_standard_errors"],wanted_qse,rtol=2e-3,atol=5e-5)
            output["assertions"]+=12+len(theta)+len(reference["standard_errors"])+1
            case["results"][layer]={"fit":result,"prediction":pred,"independent_average_nll_score":score.tolist(),"independent_score_max":score_max}
        output["cases"].append(case)
    output["confirmed_boundary_regressions"], output["regression_assertions"]=confirmed_boundary_regressions(boundaries,nullable_boundaries)
    output["assertions"]+=output["regression_assertions"]
    output["assertion_count_method"]="Scalar numerical comparisons plus explicit regression contract predicates. Nonfinite JSON is independently rejected for every decoded result."
    if args.duckdb:
        output["native_sql_queries"]=runner.queries
        output["native_sql_decoded_results"]=runner.decoded
    for flag,path in [("library",args.library),("extension",args.extension)]:
        if path:
            output[flag]={"path":str(path.resolve()),"sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"cases":len(output["cases"]),"confirmed_boundary_regressions":len(output["confirmed_boundary_regressions"]),"assertions":output["assertions"],"regression_assertions":output["regression_assertions"],"boundaries":list(boundaries),"output":str(args.output)}))

if __name__=="__main__":
    main()
