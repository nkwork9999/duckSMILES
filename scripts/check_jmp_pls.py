#!/usr/bin/env python3
"""Independent LAPACK/scikit-learn PLS references through actual C ABI and SQL.

Production uses an eigensystem in the smaller cross-covariance Gram matrix.
The regression-deflation reference is sklearn's LAPACK-SVD route, checked
against its iterative NIPALS at tol1e-26. SIMPLS uses direct LAPACK SVD and
orthogonal projections. Conditional-score uncertainty is explicitly approximate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy import stats
import sklearn
from sklearn.cross_decomposition import PLSRegression
from check_jmp_censored import Abi
from check_jmp_continuation import Runner, expression

ROOT=Path(__file__).resolve().parents[1]


def normalize(a, scale):
    origin=a[0].copy()
    delta=a-origin
    offset=delta.mean(0)
    centered=delta-offset
    sd=np.std(centered,axis=0,ddof=1) if scale else np.ones(a.shape[1])
    sd[sd==0]=1
    return centered/sd,origin,offset,sd


def reference(X,Y,Q,A,method,sx,sy):
    x,xorigin,xoffset,xsd=normalize(X,sx)
    y,yorigin,yoffset,ysd=normalize(Y,sy)
    if method==0:
        fit=PLSRegression(n_components=A,scale=False,max_iter=10000,tol=1e-26)
        fit.algorithm='svd'
        fit.fit(x,y)
        B=fit.coef_.T
        R=fit.x_rotations_
        T=fit.x_scores_
        C=fit.y_loadings_
        W=fit.x_weights_
        if X.shape[0]==12 and sx and sy:
            iterative=PLSRegression(n_components=A,scale=False,max_iter=10000,tol=1e-26).fit(x,y)
            np.testing.assert_allclose(iterative.coef_,fit.coef_,rtol=1e-10,atol=1e-11)
    else:
        S=x.T@y
        rotations=[];scores=[];loads=[];orth=[]
        for _ in range(A):
            r=np.linalg.svd(S,full_matrices=False)[0][:,0]
            t=x@r
            for oldr,oldt in zip(rotations,scores):
                z=oldt@t;r-=oldr*z;t-=oldt*z
            norm=np.linalg.norm(t);t/=norm;r/=norm
            loading=x.T@t;c=y.T@t;v=loading.copy()
            for old in orth:v-=old*(old@v)
            v/=np.linalg.norm(v);S-=np.outer(v,v@S)
            rotations.append(r);scores.append(t);loads.append(c);orth.append(v)
        R=np.array(rotations).T;T=np.array(scores).T;C=np.array(loads).T;W=R
        B=R@C.T
    z=((Q-xorigin)-xoffset)/xsd
    fitted=yorigin+(yoffset+x@B*ysd)
    predictions=yorigin+(yoffset+z@B*ysd)
    residual=((Y-yorigin)-yoffset)-x@B*ysd
    sse=np.sum(residual**2,axis=0)
    tt=np.sum(T**2,axis=0)
    ss=tt*np.sum(C**2,axis=0)
    weights=W/np.linalg.norm(W,axis=0)
    vip=np.sqrt(X.shape[1]*(weights**2@ss)/ss.sum())
    leverage=1/len(X)+np.sum((z@R)**2/tt,axis=1)
    df=len(X)-A-1
    mse=sse/df
    meanse=np.sqrt(leverage[:,None]*mse)
    indivse=np.sqrt((1+leverage[:,None])*mse)
    t=stats.t.ppf(.975,df)
    return {'coefficients':B*ysd/xsd[:,None], 'coefficients_standardized':B,
            'vip':vip,'fitted':fitted,'sse':sse,
            'predict':{'predictions':predictions,'leverage':leverage,
                       'mean_standard_error':meanse,'individual_standard_error':indivse,
                       'mean_ci_lower':predictions-t*meanse,'mean_ci_upper':predictions+t*meanse,
                       'individual_ci_lower':predictions-t*indivse,'individual_ci_upper':predictions+t*indivse}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library',type=Path,required=True)
    parser.add_argument('--duckdb')
    parser.add_argument('--extension',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    abi=Abi(args.library)
    boundaries={'exported_C_ABI':abi}
    nulls={'exported_C_ABI':abi.nullable}
    if args.duckdb:
        runner=Runner(args.duckdb,args.extension)
        boundaries['loaded_DuckDB_SQL']=lambda op,data,a:runner.one(expression(op,data,a))
        nulls['loaded_DuckDB_SQL']=lambda op,data,a:runner.query('SELECT '+expression(op,data,a)+' AS result;')[0]['result']
    checks=0;bad=0;maximum=0.;records=[]
    def check(actual,wanted,label):
        nonlocal checks,maximum
        actual=np.asarray(actual,dtype=float);wanted=np.asarray(wanted,dtype=float)
        assert actual.shape==wanted.shape,(label,actual.shape,wanted.shape)
        np.testing.assert_allclose(actual,wanted,rtol=3e-8,atol=3e-10,err_msg=label)
        checks+=actual.size
        if actual.size:maximum=max(maximum,float(np.max(np.abs(actual-wanted)/(1+np.abs(wanted)))))
    rng=np.random.default_rng(20261024)
    fixtures=[]
    for n,p,q in [(12,4,2),(8,10,3),(20,1,3)]:
        X=rng.normal(size=(n,p));Y=X@rng.normal(size=(p,q))+rng.normal(size=(n,q))*.4
        Q=rng.normal(size=(3,p));A=min(p,n-2,3)
        for method in [0,1]:
            for sx in [0,1]:
                for sy in [0,1]:fixtures.append((X,Y,Q,A,method,sx,sy))
    for index,(X,Y,Q,A,method,sx,sy) in enumerate(fixtures):
        wanted=reference(X,Y,Q,A,method,sx,sy)
        data=np.column_stack([X,Y]).tolist();params=[Y.shape[1],A,method,sx,sy]
        for layer,call in boundaries.items():
            fit=call('pls2_fit',data,params)
            assert fit['components_used']==A
            for field in ['coefficients','coefficients_standardized','vip','fitted','sse']:
                check(fit[field],wanted[field],f'{layer}/{index}/{field}')
            pred=call('pls2_predict',Q.tolist(),fit['predict_parameters'])
            for field,value in wanted['predict'].items():check(pred[field],value,f'{layer}/{index}/{field}')
            T=np.asarray(fit['x_scores']);gram=T.T@T
            check(gram-np.diag(np.diag(gram)),np.zeros_like(gram),f'{layer}/{index}/score_orthogonality')
            check(sum(v*v for v in fit['vip']),X.shape[1],f'{layer}/{index}/VIP_normalization')
            records.append({'boundary':layer,'fixture':index,'method':method,'scale_x':sx,'scale_y':sy,'coefficients':fit['coefficients'],'vip':fit['vip'],'prediction':pred})
        print('PASS regression fixture',index,flush=True)
    X=rng.normal(size=(18,4));labels=np.array([-10.,7.,100.]);groups=np.repeat(np.arange(3),6)
    X[:,0]+=groups;Y=np.eye(3)[groups];Q=np.vstack([X[:3],np.array([[8.,-8.,3.,2.]])])
    for method in [0,1]:
        for sx in [0,1]:
            for sy in [0,1]:
                wanted=reference(X,Y,Q,2,method,sx,sy)
                for layer,call in boundaries.items():
                    fit=call('plsda_fit',np.column_stack([X,labels[groups]]).tolist(),[2,method,sx,sy])
                    pred=call('plsda_predict',Q.tolist(),fit['predict_parameters'])
                    check(fit['class_labels'],labels,f'{layer}/DA/labels')
                    check(fit['coefficients'],wanted['coefficients'],f'{layer}/DA/coefficients')
                    check(fit['vip'],wanted['vip'],f'{layer}/DA/VIP')
                    check(pred['predictions'],wanted['predict']['predictions'],f'{layer}/DA/predictions')
                    check(pred['predicted_labels'],labels[np.argmax(wanted['predict']['predictions'],axis=1)],f'{layer}/DA/argmax')
                    assert 'not_calibrated_probabilities' in pred['score_semantics']
                    check(np.sum(pred['predictions'],axis=1),np.ones(len(Q)),f'{layer}/DA/indicator_sum')
    fitdata=np.column_stack([fixtures[0][0],fixtures[0][1]]).tolist()
    for layer,call in boundaries.items():
        saved=call('pls2_fit',fitdata,[2,2])['predict_parameters']
        cases=[('pls2_fit',fitdata,a) for a in [[],[0,1],[2,0],[2,2,2],[2,2,0,2],[2,2,0,1,1,1],[2,2,0,1,1,0],[2,2,0,1,1,.95,1]]]
        cases += [('pls2_predict',[[1.,2.,3.,4.]],a) for a in [[],saved[:-1],saved+[0], [99]+saved[1:],saved[:6]+[1]+saved[7:],saved[:7]+[0]+saved[8:]]]
        cases += [('pls2_predict',[[1.]],saved),('plsda_predict',[[1.,2.,3.,4.]],saved),('plsda_fit',[[0,1],[1,1]],[1])]
        for op,data,a in cases:
            assert nulls[layer](op,data,a) is None,(layer,op,data,a)
            bad+=1
        const=call('pls2_fit',[[0,2],[0,4],[0,6]],[1,1])
        assert const['components_used']==0 and const['vip'] is None
        pred=call('pls2_predict',[[100]],const['predict_parameters'])
        check(pred['predictions'],[[4]],f'{layer}/constant_x')
        check(pred['mean_standard_error'],[[np.sqrt(4/3)]],f'{layer}/intercept_SE')
        saturated=call('pls2_fit',[[0,1],[1,3]],[1,1])
        assert call('pls2_predict',[[2]],saturated['predict_parameters'])['mean_standard_error'] is None
    sources=['crates/analytics/src/lib.rs','crates/analytics/src/pls_multivariate.rs','crates/analytics/src/explore.rs','crates/analytics/src/glm.rs','crates/analytics/tests/pls_multivariate_contracts.rs','src/include/analytics_registry.inc']
    result={'schema_version':1,'status':'PASS','scope':'Centered NIPALS/SIMPLS multiresponse regression,VIP,indicator PLS-DA and conditional-score approximate prediction intervals; no installed JMP comparison','platform':platform.platform(),'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__,'boundaries':list(boundaries),'regression_fixtures_per_boundary':len(fixtures),'DA_fixtures_per_boundary':8,'checks':checks,'invalid_NULL_checks':bad,'maximum_scaled_error':maximum,'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},'library_sha256':hashlib.sha256(args.library.read_bytes()).hexdigest(),'records':records}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))

if __name__=='__main__':main()
