#!/usr/bin/env python3
"""Independent raw-likelihood ordinal oracle for the actual production C ABI.

The Rust estimator uses centered/scaled predictors, direct thresholds, analytical
scores/information and damped Newton. This test-only oracle uses original-unit
predictors, SciPy BFGS over log threshold gaps, then 80-digit numerical score
roots and Hessians of the direct multinomial logistic-CDF differences. It does
not import or translate production likelihood, gradient or quantile helpers.
"""
import argparse
import ctypes
import datetime
import hashlib
import json
import math
from pathlib import Path
import platform
import traceback

import mpmath as mp
import numpy as np
import scipy
from scipy import optimize, special, stats

ROOT = Path(__file__).resolve().parents[1]
CHECKS = 0
MAX_SCALED_ERROR = 0.0
REFERENCES = [
    "https://www.jmp.com/support/help/en/19.2/jmp/ordinal-responses.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/specify-a-distribution.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/parameter-estimates-for-original-predictors.shtml",
    "https://docs.scipy.org/doc/scipy/reference/optimize.minimize-bfgs.html",
]


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def close(actual, expected, label, atol=2e-7, rtol=2e-7):
    global MAX_SCALED_ERROR
    a = np.asarray(actual, dtype=float)
    e = np.asarray(expected, dtype=float)
    require(a.shape == e.shape, (label, "shape", a.shape, e.shape))
    require(np.isfinite(a).all(), (label, "nonfinite", actual))
    error = np.abs(a-e)
    scaled = np.max(error/(1+np.abs(e))) if error.size else 0.
    MAX_SCALED_ERROR = max(MAX_SCALED_ERROR, float(scaled))
    require(np.all(error <= atol+rtol*np.abs(e)), (label, actual, np.asarray(expected).tolist(), error.tolist()))


def strict_json(raw):
    def invalid(value):
        raise AssertionError(("non-JSON numeric constant", value))
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, ("duplicate JSON key", key))
            result[key] = value
        return result
    result = json.loads(raw, parse_constant=invalid, object_pairs_hook=unique)
    def finite(value):
        if isinstance(value, (float, int)) and not isinstance(value, bool):
            require(math.isfinite(value), ("nonfinite JSON value", value))
        elif isinstance(value, list):
            for item in value:
                finite(item)
        elif isinstance(value, dict):
            for item in value.values():
                finite(item)
    finite(result)
    return result


class Abi:
    def __init__(self, path):
        self.library = ctypes.CDLL(str(path))
        double_p = ctypes.POINTER(ctypes.c_double)
        self.library.ds_analytics_run.argtypes = [
            ctypes.c_char_p, ctypes.c_size_t, double_p, ctypes.c_size_t,
            ctypes.c_size_t, double_p, ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_size_t)]
        self.library.ds_analytics_run.restype = ctypes.c_void_p
        self.library.ds_analytics_free.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        self.library.ds_analytics_free.restype = None
        self.calls = 0

    def __call__(self, operation, rows, args=()):
        data = np.ascontiguousarray(rows, dtype=np.float64)
        parameters = np.ascontiguousarray(args, dtype=np.float64)
        require(data.ndim == 2, "oracle rectangular data")
        encoded = operation.encode()
        length = ctypes.c_size_t(731)
        ptr = self.library.ds_analytics_run(
            encoded, len(encoded), data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            *data.shape, parameters.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            len(parameters), ctypes.byref(length))
        self.calls += 1
        if not ptr:
            require(length.value == 0, "public NULL must clear output length")
            return None
        try:
            return strict_json(ctypes.string_at(ptr, length.value))
        finally:
            self.library.ds_analytics_free(ptr, length.value)


def probability(theta, features, levels):
    """Direct high-precision CDF subtraction in original model coordinates."""
    m = levels-1
    eta = sum(a*b for a, b in zip(theta[m:], features))
    cumulative = [1/(1+mp.exp(-a-eta)) for a in theta[:m]]
    return [cumulative[0]] + [b-a for a, b in zip(cumulative, cumulative[1:])] + [1-cumulative[-1]]


def likelihood(rows, labels):
    p = len(rows[0])-2
    data = [([mp.mpf(x) for x in row[:p]], labels.index(row[p]), mp.mpf(row[p+1]))
            for row in rows if row[-1] > 0]
    def ll(*theta):
        return sum(weight*mp.log(probability(theta, x, len(labels))[y])
                   for x, y, weight in data)
    return ll


def derivative(function, theta, i, j=None):
    order = [0]*len(theta)
    order[i] += 1
    if j is not None:
        order[j] += 1
    return mp.diff(function, tuple(theta), tuple(order))


def direct_reference(name, rows, known=None):
    labels = sorted({r[-2] for r in rows if r[-1] > 0})
    p, m = len(rows[0])-2, len(labels)-1
    ll = likelihood(rows, labels)
    if known is None:
        counts = np.array([sum(r[-1] for r in rows if r[-2] == label) for label in labels])
        cuts = special.logit(np.cumsum(counts)[:-1]/counts.sum())
        initial = np.r_[cuts[0], np.log(np.diff(cuts)), np.zeros(p)]
        def unpack(t):
            return np.r_[np.r_[t[0], t[0]+np.cumsum(np.exp(t[1:m]))], t[m:]]
        def objective(t):
            raw = unpack(t)
            value = 0.
            for row in rows:
                c = special.expit(raw[:m]+np.dot(raw[m:], row[:p]))
                probabilities = np.diff(np.r_[0., c, 1.])
                q = probabilities[labels.index(row[p])]
                if q <= 0:
                    return float("inf")
                value -= row[-1]*np.log(q)
            return value
        starts = [initial, initial+np.r_[np.zeros(m), np.full(p, .25)],
                  initial+np.r_[np.zeros(m), np.full(p, -.25)]]
        fits = [optimize.minimize(objective, start, method="BFGS", jac="3-point",
                                 options={"gtol":1e-9,"maxiter":1000}) for start in starts]
        best = min(fits, key=lambda r:r.fun)
        first = [mp.mpf(x) for x in unpack(best.x)]
        def scores(*theta):
            return tuple(derivative(ll, theta, i) for i in range(len(theta)))
        theta = list(mp.findroot(scores, tuple(first), tol=mp.mpf("1e-60"), maxsteps=30))
        optimizer = [{"method":"BFGS; numerical3point derivative; log threshold gaps",
                      "success":bool(r.success),"objective":float(r.fun),"message":str(r.message)} for r in fits]
    else:
        theta = list(known)
        optimizer = "independent exact parameter fixture; no optimizer"
    score = [derivative(ll, theta, i) for i in range(len(theta))]
    require(max(abs(s) for s in score) < mp.mpf("1e-50"), (name,"reference score not zero",[str(s) for s in score]))
    information = mp.matrix([[-derivative(ll,theta,i,j) for j in range(len(theta))]
                             for i in range(len(theta))])
    covariance = information**-1
    require(min(np.linalg.eigvalsh(np.array(information.tolist(),dtype=float))) > 0,
            (name,"reference information not PD"))
    transform = mp.eye(len(theta))
    for j in range(1,m):
        transform[j,j-1] = -1
    genreg_cov = transform*covariance*transform.T
    genreg = [theta[0]] + [theta[j]-theta[j-1] for j in range(1,m)]
    return {"name":name,"rows":rows,"labels":labels,"p":p,"K":len(labels),
            "theta":theta,"covariance":covariance,"genreg":genreg,
            "genreg_covariance":genreg_cov,"log_likelihood":ll(*theta),
            "optimizer":optimizer,"score_max":max(abs(s) for s in score)}


def fixtures():
    mp.mp.dps = 80
    intercept = direct_reference("three_class_intercept_N100", [[10.,40.],[20.,35.],[30.,25.]],
                                 [mp.log(mp.mpf(2)/3),mp.log(3)])
    exact_cov = [[mp.mpf(1)/24,mp.mpf(1)/45],[mp.mpf(1)/45,mp.mpf(4)/75]]
    close(intercept["covariance"].tolist(),exact_cov,"analytical multinomial delta covariance",1e-15,1e-15)
    grouped = [[x,label,float(n)] for x,counts in [(-1.,[2,8,10]),(0.,[5,10,5]),(1.,[10,8,2])]
               for label,n in zip([10.,20.,30.],counts)]
    symmetric = direct_reference("exact_common_slope_sign",grouped,[-mp.log(3),mp.log(3),mp.log(3)])
    shifted = [[10.+3.*r[0],r[1],r[2]] for r in grouped]
    beta = mp.log(3)/3
    translated = direct_reference("translated_scaled_predictor",shifted,[-mp.log(3)-10*beta,mp.log(3)-10*beta,beta])
    cases = [(-1.,-1.,[2,5,9,14]),(-1.,1.,[5,8,10,7]),
             (0.,-1.,[5,8,10,7]),(0.,1.,[8,10,8,4]),
             (1.,-1.,[10,9,7,4]),(1.,1.,[15,8,5,2])]
    multiple = [[x,z,label,float(n)] for x,z,counts in cases
                for label,n in zip([-10.,0.,4.,90.],counts)]
    general = direct_reference("two_predictors_four_levels",multiple)
    require(all(b>0 for b in general["theta"][3:]),"both cumulative slopes positive in independent fixture")
    return [intercept,symmetric,translated,general]


def mp_lists(matrix):
    return [[float(v) for v in row] for row in matrix.tolist()]


def prediction_reference(theta, covariance, query, labels, confidence):
    p = len(query)-1
    probabilities = probability(theta,[mp.mpf(x) for x in query[:p]],len(labels))
    cumulative = [sum(probabilities[:j+1]) for j in range(len(labels))]
    q = mp.mpf(query[-1])
    chosen = next(j for j,v in enumerate(cumulative) if v>=q)
    mode = max(range(len(labels)),key=lambda i:probabilities[i])
    # Only queries away from ties are used for binary64 fit comparisons.
    require(min(abs(v-q) for v in cumulative[:-1]) > mp.mpf("1e-6"),"reference query away from CDF ties")
    errors=[]
    for j in range(len(labels)):
        def f(*t):return probability(t,[mp.mpf(x) for x in query[:p]],len(labels))[j]
        g=mp.matrix([derivative(f,theta,i) for i in range(len(theta))])
        errors.append(mp.sqrt((g.T*covariance*g)[0]))
    z = stats.norm.ppf((1+confidence)/2)
    return {"probabilities":[float(v) for v in probabilities],
            "cumulative_probabilities":[float(v) for v in cumulative],
            "most_likely":labels[mode],"quantiles":labels[chosen],
            "linear_predictor":float(sum(a*mp.mpf(b) for a,b in zip(theta[len(labels)-1:],query[:p]))),
            "probability_standard_errors":[float(v) for v in errors],
            "probability_confidence_intervals":[[float(v-z*s),float(v+z*s)] for v,s in zip(probabilities,errors)]}


def check_fit(abi, ref):
    fit=abi("genreg_ordinal_fit",ref["rows"],[1000.,1e-10,.95])
    require(fit is not None,(ref["name"],"fit NULL"))
    require(fit["converged"] and fit["inference_valid"],(ref["name"],"fit inference invalid"))
    close(fit["cutpoints"],ref["theta"][:ref["K"]-1],"raw cutpoints")
    close(fit["coefficients"],ref["theta"][ref["K"]-1:],"original slopes")
    close(fit["covariance"],mp_lists(ref["covariance"]),"observed information covariance")
    close(fit["genreg_intercepts"],ref["genreg"],"first cut and differences")
    close(fit["genreg_covariance"],mp_lists(ref["genreg_covariance"]),"display covariance")
    close(fit["log_likelihood"],ref["log_likelihood"],"full multinomial ll")
    d=len(ref["theta"]);n=sum(r[-1] for r in ref["rows"])
    close(fit["aic"],2*d-2*ref["log_likelihood"],"aic")
    close(fit["bic"],d*math.log(n)-2*ref["log_likelihood"],"bic effective N")
    require(fit["covariance_method"]=="observed_information" and not fit["jmp_sandwich_equivalence_verified"],"named uncertainty residual")
    z=stats.norm.ppf(.975)
    se=np.sqrt(np.diag(mp_lists(ref["covariance"])))
    close(fit["standard_errors"],se,"SE")
    theta=np.array([float(v) for v in ref["theta"]])
    close(fit["confidence_intervals"],np.column_stack([theta-z*se,theta+z*se]),"Wald CI")
    carrier=fit["predict_parameters"]
    require(len(carrier)==5+3*ref["p"]+ref["K"]+d+d*d,"exact carrier length")
    require(len(carrier)<=1061,"public carrier cap")
    if ref["p"]==0:queries=[[.13],[.63],[.91]]
    elif ref["name"]=="translated_scaled_predictor":queries=[[7.,.13],[10.,.63],[13.,.91]]
    elif ref["p"]==1:queries=[[-1.,.13],[0.,.63],[1.,.91]]
    else:queries=[[-.7,.3,.13],[.2,-.8,.63],[1.4,.2,.91]]
    prediction=abi("genreg_ordinal_predict",queries,carrier)
    require(prediction is not None,"saved prediction NULL")
    references=[prediction_reference(ref["theta"],ref["covariance"],q,ref["labels"],.95) for q in queries]
    for key in references[0]:close(prediction[key],[r[key] for r in references],"replay "+key)
    single=abi("genreg_ordinal_predict",queries[:1],carrier)
    close(single["probabilities"][0],prediction["probabilities"][0],"batch invariance",0,0)
    return {"name":ref["name"],"fit":fit,"queries":queries,"prediction":prediction,
            "prediction_reference":references}


def invalid_cases(abi, carrier):
    cases=[]
    def rejected(name,op,rows,args):
        require(abi(op,rows,args) is None,(name,"malformed request accepted"))
        cases.append(name)
    predict="genreg_ordinal_predict"
    for name,change in [
        ("carrier truncated",lambda a:a[:-1]),("carrier trailing value",lambda a:a+[0.]),
        ("bad version",lambda a:[2.]+a[1:]),("too many predictors",lambda a:a[:1]+[17.]+a[2:]),
        ("too many levels",lambda a:a[:2]+[17.]+a[3:]),("bad confidence",lambda a:a[:3]+[1.]+a[4:]),
        ("nonboolean covariance",lambda a:a[:4]+[2.]+a[5:])]:
        rejected(name,predict,[[.63]],change(list(carrier)))
    # The intercept-only carrier has p0,K3, labels at5..8,theta8..10,cov10..14.
    for name,index,value in [("unordered labels",6,10.),("unordered cuts",9,-1.),
                             ("negative variance",10,-1.),("asymmetric covariance",11,.1)]:
        a=list(carrier);a[index]=value;rejected(name,predict,[[.63]],a)
    a=list(carrier);a[4]=0.;rejected("disabled covariance nonzero",predict,[[.63]],a)
    for q in [0.,1.,-.1]:rejected("invalid quantile "+str(q),predict,[[q]],carrier)
    fit="genreg_ordinal_fit"
    for name,rows in [
        ("fractional frequency",[[10.,.5],[20.,2.]]),
        ("single effective level",[[10.,4.],[20.,0.]]),
        ("seventeen levels",[[float(i),1.] for i in range(17)]),
        ("per row frequency cap",[[10.,1000001.],[20.,2.]]),
        ("total frequency cap",[[float(i%2),1000000.] for i in range(1001)]),
        ("constant predictor",[[1.,10.,2.],[1.,20.,2.]]),
        ("finite label magnitude cap",[[1e151,1.],[0.,1.]])]:rejected(name,fit,rows,[])
    # No covariance and 16 outcomes; reject the cell bound before row evaluation.
    a=[1.,0.,16.,.95,0.]+list(map(float,range(16)))+list(map(float,range(-7,8)))+[0.]*225
    rejected("prediction probability cells cap",predict,[[.63]]*12501,a)
    return cases


def serializable_reference(ref):
    return {"name":ref["name"],"data":ref["rows"],"ordered_labels":ref["labels"],
            "raw_parameters_80digit":[mp.nstr(v,75) for v in ref["theta"]],
            "raw_covariance_80digit":[[mp.nstr(v,75) for v in row] for row in ref["covariance"].tolist()],
            "genreg_parameters_80digit":[mp.nstr(v,75) for v in ref["genreg"]],
            "genreg_covariance_80digit":[[mp.nstr(v,75) for v in row] for row in ref["genreg_covariance"].tolist()],
            "log_likelihood_80digit":mp.nstr(ref["log_likelihood"],75),
            "score_max":str(ref["score_max"]),"optimizer":ref["optimizer"]}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--library",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--reference-only",action="store_true")
    args=parser.parse_args()
    require(args.reference_only or args.library is not None,"library required for runtime oracle")
    evidence={"schema":"jmp-ordinal-independent-reference-v1","status":"RUNNING",
              "time_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
              "mpmath":mp.__version__,"precision_digits":80,"primary_references":REFERENCES,
              "source_sha256":hashlib.sha256((ROOT/"crates/analytics/src/genreg_ordinal.rs").read_bytes()).hexdigest(),
              "harness_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "comparison":"raw original-unit direct multinomial likelihood; observed information only, not JMP sandwich equivalence",
              "installed_JMP_comparison":False,"whole_JMP_family_parity":False}
    try:
        refs=fixtures()
        evidence["references"]=[serializable_reference(r) for r in refs]
        if not args.reference_only:
            abi=Abi(args.library)
            evidence["library_sha256"]=hashlib.sha256(args.library.read_bytes()).hexdigest()
            evidence["public_c_abi"]=[check_fit(abi,r) for r in refs]
            evidence["invalid_NULL_cases"]=invalid_cases(abi,evidence["public_c_abi"][0]["fit"]["predict_parameters"])
            evidence["calls"]=abi.calls
        evidence["status"]="PASS_REFERENCE_ONLY" if args.reference_only else "PASS"
    except Exception:
        evidence["status"]="FAIL"
        evidence["exception"]=traceback.format_exc()
        raise
    finally:
        evidence["checks"]=CHECKS
        evidence["maximum_scaled_error"]=MAX_SCALED_ERROR
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(evidence,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:evidence[k] for k in ["status","checks","maximum_scaled_error"]}))

if __name__=="__main__":main()
