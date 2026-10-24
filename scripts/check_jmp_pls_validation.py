#!/usr/bin/env python3
"""Independent per-split sklearn/LAPACK PLS selection and mpmath underflow oracle."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys

import mpmath as mp
import numpy as np
import scipy
import sklearn
from check_jmp_pls import reference,normalize
from check_jmp_censored import Abi
from check_jmp_continuation import Runner,expression

ROOT=Path(__file__).resolve().parents[1]


def oracle(X,Y,ids,A,method,sx,sy,mode,da=False):
    ref=np.arange(len(X)) if mode==0 else np.flatnonzero(ids==0)
    valid=np.arange(len(X)) if mode==0 else np.flatnonzero(ids==1)
    test=np.array([],int) if mode==0 else np.flatnonzero(ids==2)
    groups=[np.flatnonzero(ids==f) for f in sorted(set(ids))] if mode==0 else [valid]
    variance=Y[ref].var(0,ddof=1)
    path=[];press=[];root=[];standardized=[];q2=[];transforms=[]
    for a in range(A+1):
        pred=np.zeros((len(valid),Y.shape[1]));sse=np.zeros(Y.shape[1]);std=np.zeros(Y.shape[1])
        for index,v in enumerate(groups):
            tr=np.setdiff1d(ref,v) if mode==0 else ref
            values=np.tile(Y[tr].mean(0),(len(v),1)) if a==0 else reference(X[tr],Y[tr],X[v],a,method,sx,sy)['predict']['predictions']
            pred[np.searchsorted(valid,v)]=values
            errors=(Y[v]-values)**2;sse+=errors.sum(0);std+=errors.mean(0)/variance
            if a==0:
                _,xo,xm,xs=normalize(X[tr],sx);_,yo,ym,ys=normalize(Y[tr],sy)
                if da:yo=np.zeros(Y.shape[1]);ym=Y[tr].mean(0);ys=np.std(Y[tr],axis=0,ddof=1) if sy else np.ones(Y.shape[1])
                transforms.append({'group':int(ids[v[0]]) if mode==0 else 1,'training_count':len(tr),'validation_count':len(v),'feature_origins':xo,'feature_mean_offsets':xm,'feature_scales':xs,'response_origins':yo,'response_mean_offsets':ym,'response_scales':ys})
        std/=max(1,len(groups)-1)
        path.append(pred);press.append(sse);standardized.append(std);root.append(np.sqrt(std.mean()))
        sst=((Y[valid]-Y[valid].mean(0))**2).sum(0);q2.append(1-sse/sst)
    selected=int(np.argmin(root));final=np.tile(Y[ref].mean(0),(len(X),1)) if selected==0 else reference(X[ref],Y[ref],X,selected,method,sx,sy)['predict']['predictions']
    return {'predictions_by_components':np.array(path),'raw_press_per_response':press,'normalized_press_per_response':standardized,'root_mean_press':root,'pooled_predictive_r2_per_response':q2,'reference_response_variances':variance,'selected_components':selected,'validation_indices':valid,'test_indices':test,'test_predictions':final[test] if len(test) else np.empty((0,)),'final_predictions':final,'split_training_transforms':transforms}


def mp_ols(x,y,query):
    xm=sum(x)/len(x);ym=sum(y)/len(y)
    slope=sum((a-xm)*(b-ym) for a,b in zip(x,y))/sum((a-xm)**2 for a in x)
    return [ym+slope*(a-xm) for a in query]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library',type=Path,required=True);parser.add_argument('--duckdb');parser.add_argument('--extension',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();abi=Abi(args.library)
    calls={'exported_C_ABI':abi};nulls={'exported_C_ABI':abi.nullable}
    if args.duckdb:
        runner=Runner(args.duckdb,args.extension)
        calls['loaded_DuckDB_SQL']=lambda op,d,a:runner.one(expression(op,d,a))
        nulls['loaded_DuckDB_SQL']=lambda op,d,a:runner.query('SELECT '+expression(op,d,a)+' AS result;')[0]['result']
    checks=0;bad=0;maximum=0.;records=[]
    def compare(actual,wanted,label,rtol=3e-8,atol=3e-10):
        nonlocal checks,maximum
        a=np.asarray(actual,dtype=float);b=np.asarray(wanted,dtype=float)
        assert a.shape==b.shape,(label,a.shape,b.shape)
        np.testing.assert_allclose(a,b,rtol=rtol,atol=atol,err_msg=label)
        checks+=a.size
        if a.size:maximum=max(maximum,float(np.max(np.abs(a-b)/(1+np.abs(b)))))
    rng=np.random.default_rng(20261025);fixtures=[]
    for mode,n,p,q,A in [(0,12,4,2,3),(0,8,3,1,2),(1,15,3,2,2)]:
        X=rng.normal(size=(n,p))*np.arange(1,p+1)+np.arange(p)*3
        Y=X@rng.normal(size=(p,q))+rng.normal(size=(n,q))*.5
        ids=(np.array([0]*2+[1]*4+[2]*6) if n==12 else np.arange(n)) if mode==0 else np.array([0]*8+[1]*4+[2]*3)
        for method in [0,1]:
            for sx,sy in [(0,0),(1,1)]:fixtures.append((False,X,Y,ids,A,method,sx,sy,mode))
    for mode in [0,1]:
        n=18;X=rng.normal(size=(n,4));groups=np.tile(np.arange(3),6);X[:,0]+=groups;Y=np.eye(3)[groups]
        ids=np.arange(n)%2 if mode==0 else np.array([0]*9+[1]*6+[2]*3)
        for method in [0,1]:fixtures.append((True,X,Y,ids,2,method,1,1,mode))
    labels=np.array([-10.,7.,100.])
    for index,(da,X,Y,ids,A,method,sx,sy,mode) in enumerate(fixtures):
        wanted=oracle(X,Y,ids,A,method,sx,sy,mode,da)
        data=np.column_stack([X,labels[np.argmax(Y,axis=1)] if da else Y,ids]).tolist()
        params=([A,method,sx,sy,mode] if da else [Y.shape[1],A,method,sx,sy,mode]);op='plsda_validate' if da else 'pls2_validate'
        for boundary,call in calls.items():
            actual=call(op,data,params)
            assert actual['selected_components']==wanted['selected_components']
            for field in ['predictions_by_components','raw_press_per_response','normalized_press_per_response','root_mean_press','pooled_predictive_r2_per_response','reference_response_variances','validation_indices','test_indices','test_predictions']:
                compare(actual[field],wanted[field],f'{boundary}/{index}/{field}')
            for a,b in zip(actual['split_training_transforms'],wanted['split_training_transforms']):
                for field,value in b.items():compare(a[field],value,f'{boundary}/{index}/provenance/{field}')
            saved=actual['predict_parameters'];assert saved[0]==2
            predicted=call('plsda_predict' if da else 'pls2_predict',X.tolist(),saved)
            compare(predicted['predictions'],wanted['final_predictions'],f'{boundary}/{index}/selected_final_model')
            if da:compare(actual['selected_validation_labels'],labels[np.argmax(wanted['predictions_by_components'][wanted['selected_components']],axis=1)],f'{boundary}/{index}/DAlabels')
            if mode==1:
                changed=np.array(data);changed[ids==2,0]+=1e6;changed[ids==2,X.shape[1]:-1]+=0 if da else 1e6
                alternate=call(op,changed.tolist(),params)
                compare(alternate['root_mean_press'],actual['root_mean_press'],f'{boundary}/{index}/test_ignored',0,0)
                compare(alternate['predict_parameters'],saved,f'{boundary}/{index}/test_not_trained',0,0)
            records.append({'boundary':boundary,'fixture':index,'DA':da,'mode':mode,'method':method,'selected_components':actual['selected_components'],'root_mean_press':actual['root_mean_press'],'predictions_by_components':actual['predictions_by_components']})
        print('PASS validation fixture',index,flush=True)
    mp.mp.dps=90;noise=np.array([0,1,-1,2,-2,1.]);X=np.arange(6,dtype=float);Y=(2*X+noise*1e-8)*1e-155;ids=np.arange(6)%2
    x=list(map(mp.mpf,map(float,X)));y=list(map(mp.mpf,map(float,Y)));ym=sum(y)/6;variance=sum((v-ym)**2 for v in y)/5;stdloss=mp.mpf(0)
    for f in [0,1]:
        tr=np.flatnonzero(ids!=f);v=np.flatnonzero(ids==f);pred=mp_ols([x[i] for i in tr],[y[i] for i in tr],[x[i] for i in v]);stdloss+=sum((y[i]-p)**2 for i,p in zip(v,pred))/len(v)/variance
    mp_root=float(mp.sqrt(stdloss));fitted=mp_ols(x,y,x);rss=sum((a-b)**2 for a,b in zip(y,fitted));mp_se=float(mp.sqrt(rss/4/6));edges={}
    for boundary,call in calls.items():
        result=call('pls2_validate',np.column_stack([X,Y,ids]).tolist(),[1,1]);compare(result['root_mean_press'][1],mp_root,f'{boundary}/mp90_tiny_root',1e-6,0)
        fit=call('pls2_fit',np.column_stack([X,Y]).tolist(),[1,1]);pred=call('pls2_predict',[[2.5]],fit['predict_parameters']);compare(pred['mean_standard_error'][0][0],mp_se,f'{boundary}/mp90_tiny_SE',1e-6,0)
        legacy=json.loads((ROOT/'tasks/jmp-pls-api.json').read_text())['operations'][1]['fixture']['args'];assert legacy[0]==1;old=call('pls2_predict',[[1]],legacy);compare(old['predictions'],[[3.6,.9]],f'{boundary}/legacy_v1_prediction');compare(old['mean_standard_error'],[[np.sqrt(.08),np.sqrt(.07)]],f'{boundary}/legacy_v1_SE')
        edges[boundary]={'mp90_root':mp_root,'actual_root':result['root_mean_press'][1],'mp90_standard_error':mp_se,'actual_standard_error':pred['mean_standard_error'][0][0],'carrier_version':fit['predict_parameters'][0]}
        base=[[0,1,0],[1,3,1],[2,5,0],[3,7,1]]
        cases=[('pls2_validate',base,a) for a in [[],[0,1],[1,0],[1,1,2],[1,1,0,2],[1,1,0,1,1,2],[1,1,0,1,1,0,1]]]
        cases += [('pls2_validate',[[0,1,0],[1,3,0],[2,5,0]],[1,1]),('pls2_validate',[[0,1,0],[1,3,1]],[1,1]),('pls2_validate',[[0,1,0],[1,1,1],[2,1,0],[3,1,1]],[1,1]),('pls2_validate',[[0,1,.5],[1,3,1],[2,5,0],[3,7,1]],[1,1]),('pls2_validate',[[0,1,0],[1,3,0],[2,5,0]],[1,1,0,1,1,1]),('plsda_validate',[[0,10,0],[1,10,0],[2,30,1],[3,30,1]],[1]),('plsda_validate',[[0,10,0],[1,10,0],[2,30,1],[3,30,2]],[1,0,1,1,1])]
        for op,data,params in cases:
            assert nulls[boundary](op,data,params) is None,(boundary,op,params)
            bad+=1
    sources=['crates/analytics/src/lib.rs','crates/analytics/src/pls_multivariate.rs','crates/analytics/tests/pls_multivariate_contracts.rs','crates/analytics/tests/pls_validation_contracts.rs','src/include/analytics_registry.inc']
    result={'schema_version':1,'status':'PASS','scope':'Supplied folds/roles,0..A component selection,training-only transforms,held-out/test isolation,stable normalized PRESS and v2 SD scoring with legacyv1 compatibility; not full JMP equivalence','platform':platform.platform(),'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__,'mpmath':mp.__version__,'boundaries':list(calls),'fixtures_per_boundary':len(fixtures),'checks':checks,'invalid_NULL_checks':bad,'maximum_scaled_error':maximum,'edge_cases':edges,'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},'library_sha256':hashlib.sha256(args.library.read_bytes()).hexdigest(),'records':records}
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))

if __name__=='__main__':main()
