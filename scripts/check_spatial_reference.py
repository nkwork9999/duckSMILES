#!/usr/bin/env python3
"""Independent dense SciPy/NumPy numerical reference for bounded spatial methods."""
import argparse,datetime,hashlib,json,math,platform,subprocess
from pathlib import Path
import numpy as np
import scipy
from scipy import linalg
from scipy.spatial.distance import cdist,pdist,squareform
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--wrapper',required=True);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
proc=subprocess.Popen([a.wrapper],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
rng=np.random.default_rng(20261025);calls=0;fields=0;maxerr=0.;groups={}
def run(op,data,args=()):
 global calls
 z=np.asarray(data);assert z.ndim==2
 proc.stdin.write(' '.join(map(str,[op,*z.shape,len(args),*z.ravel(),*args]))+'\n');proc.stdin.flush();s=proc.stdout.readline();assert s,op;calls+=1
 return json.loads(s,parse_constant=lambda v:(_ for _ in ()).throw(AssertionError(v)))
def equal(x,y,label,tol=2e-7):
 global fields,maxerr
 x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float);assert x.shape==y.shape,(label,x.shape,y.shape)
 err=np.abs(x-y)/(1+np.abs(y));error=float(np.max(err,initial=0.));assert error<tol,(label,error,x,y)
 maxerr=max(maxerr,error);fields+=x.size;groups[label]=groups.get(label,0)+x.size
def corr(h,kind,rho):
 if kind==0:
  z=h/rho;return np.where(z<1,1-1.5*z+.5*z**3,0.)
 if kind==1:return np.exp(-h/rho)
 if kind==2:return np.exp(-(h/rho)**2)
 return np.power(rho,h)
for iteration in range(120):
 kind=iteration%4;d=1+iteration%3;n=8+iteration%17;pfix=iteration%4
 loc=rng.normal(size=(n,d));features=rng.normal(size=(n,pfix));x=np.column_stack([np.ones(n),features]);y=x@rng.normal(size=pfix+1)+rng.normal(size=n)
 sill=10**rng.uniform(-3,3);nugget=sill*10**rng.uniform(-1.5,.5);rho=rng.uniform(.1,.9)if kind==3 else 10**rng.uniform(-.5,.7)
 h=squareform(pdist(loc));K=sill*corr(h,kind,rho)+nugget*np.eye(n)
 cov=run('spatial_covariance',loc,[kind,sill,nugget,rho]);assert cov is not None
 equal(cov['covariance'],K,'covariance');equal(cov['semivariance'],sill+nugget-K,'record_semivariance')
 data=np.column_stack([loc,features,y]);fit=run('spatial_gls',data,[d,pfix,kind,sill,nugget,rho]);assert fit is not None
 L=linalg.cholesky(K,lower=True);wx=linalg.solve_triangular(L,x,lower=True);wy=linalg.solve_triangular(L,y,lower=True)
 Q,R=linalg.qr(wx,mode='economic');b=linalg.solve_triangular(R,Q.T@wy);ri=linalg.solve_triangular(R,np.eye(pfix+1));bcov=ri@ri.T;res=y-x@b;alpha=linalg.cho_solve((L,True),res)
 ll=-.5*(n*np.log(2*np.pi)+2*np.log(np.diag(L)).sum()+res@alpha)
 equal([fit['intercept'],*fit['coefficients']],b,'GLS_beta');equal(fit['coefficient_covariance'],bcov,'known_covariance_beta_covariance');equal(fit['marginal_fitted'],x@b,'GLS_marginal');equal(fit['conditional_fitted'],x@b+(K-nugget*np.eye(n))@alpha,'GLS_conditional');equal(fit['log_likelihood'],ll,'known_covariance_loglikelihood')
 qloc=np.vstack([loc[:2],rng.normal(size=(4,d))]);qfeat=rng.normal(size=(6,pfix));qx=np.column_stack([np.ones(6),qfeat]);cross=sill*corr(cdist(qloc,loc),kind,rho);solve=linalg.cho_solve((L,True),cross.T);delta=qx-solve.T@x
 pred=run('spatial_predict',np.column_stack([qloc,qfeat]),fit['predict_parameters']);assert pred is not None
 var=sill-np.sum(cross*solve.T,axis=1)+np.einsum('ij,jk,ik->i',delta,bcov,delta)
 equal(pred['predictions'],qx@b+cross@alpha,'universal_kriging_mean');equal(pred['latent_variances'],var,'universal_kriging_latent_variance');equal(pred['observation_variances'],var+nugget,'independent_query_nugget_variance')
 if iteration%3==0:
  sy=10**rng.uniform(-25,25);shift=2*sy;changed=data.copy();changed[:,-1]=sy*y+shift
  f2=run('spatial_gls',changed,[d,pfix,kind,sill,nugget,rho]);assert f2 is not None
  pred2=run('spatial_predict',np.column_stack([qloc,qfeat]),f2['predict_parameters']);assert pred2 is not None
  equal(np.asarray(pred2['predictions'])/sy,pred['predictions']+np.ones(6)*2,'response_unit_equivariance');equal(pred2['latent_variances'],var,'known_covariance_not_scaled_with_response')
 if iteration%4==0:
  fx=10**rng.uniform(-20,20,size=pfix);changed=data.copy();changed[:,d:d+pfix]*=fx
  f2=run('spatial_gls',changed,[d,pfix,kind,sill,nugget,rho]);assert f2 is not None
  pred2=run('spatial_predict',np.column_stack([qloc,qfeat*fx]),f2['predict_parameters']);assert pred2 is not None
  equal(pred2['predictions'],pred['predictions'],'feature_unit_equivariance');equal(pred2['latent_variances'],var,'feature_unit_variance_equivariance')
for iteration in range(80):
 d=1+iteration%3;n=3+iteration%15;loc=rng.integers(-3,4,size=(n,d)).astype(float);y=rng.normal(size=n);bins=1+iteration%10
 distances=pdist(loc);diff=pdist(y[:,None]);lo=float(distances.min());hi=float(distances.max());nb=1 if lo==hi else bins
 edges=np.linspace(lo,hi,nb+1);idx=np.minimum(np.searchsorted(edges,distances,side='right')-1,nb-1)
 out=run('spatial_variogram',np.column_stack([loc,y]),[bins]);assert out is not None
 counts=np.bincount(idx,minlength=nb);equal(out['pair_counts'],counts,'empirical_pair_counts');equal(out['edges'],edges,'empirical_edges')
 for k,c in enumerate(counts):
  if c:
   equal(out['mean_distances'][k],distances[idx==k].mean(),'empirical_mean_distance');equal(out['semivariances'][k],.5*np.mean(diff[idx==k]**2),'Matheron_semivariance')
  else:assert out['semivariances'][k]is None and out['mean_distances'][k]is None
# Duplicate records have independent measurement errors; process cross covariance excludes the nugget.
for kind in range(4):
 rho=.5 if kind==3 else 2.;out=run('spatial_covariance',[[0],[0],[1]],[kind,2.,3.,rho]);equal(out['covariance'][0],[5.,2.,2*corr(np.array([1.]),kind,rho)[0]],'collocated_records')
for data,args in [([[0,1],[0,2],[1,3]],[1,0,1,1,0,1]),([[0,1,1],[1,1,2],[2,1,3]],[1,1,1,1,1,1])]:assert run('spatial_gls',data,args)is None
# Independent analytic limits avoid subtracting rounded covariance entries.
for kind,d,rho,want in [(0,1e-20,1,1.5e-20),(1,1e-20,1,1e-20),(2,1e-10,1,1e-20),(3,1e-20,.5,math.log(2)*1e-20)]:
 out=run('spatial_covariance',[[0],[d]],[kind,1,0,rho]);equal(out['semivariance'][0][1]/want,1.,'small_distance_semivariance')
out=run('spatial_covariance',[[0],[0]],[1,1e16,1,1]);equal(out['semivariance'][0][1],1.,'independent_tiny_nugget')
for kind,d,rho,want in [(0,5e-324,1e100,1e-323),(1,5e-324,1e100,5e-324),(2,1e-200,1,1e-300),(3,5e-324,.5,5e-324*1e100*math.log(2))]:
 out=run('spatial_covariance',[[0],[d]],[kind,1e100,0,rho]);equal(out['semivariance'][0][1]/want,1.,'physical_semivariance_underflow')
for sill in [1.,1e100]:
 fit=run('spatial_gls',[[0,1],[1,2],[2,3]],[1,0,1,sill,0,1]);out=run('spatial_predict',[[0],[1],[2]],fit['predict_parameters']);equal(out['predictions'],[1,2,3],'noiseless_interpolation');assert out['latent_variances']==[0]*3 and out['observation_variances']==[0]*3
out=run('spatial_gls',[[0,-1e-200,-1e-200],[1,0,0],[2,1e-200,1e-200]],[1,1,0,1e-100,0,1]);equal(out['coefficients'],[1.],'finite_physical_beta');equal(out['coefficient_covariance'][1][1]/5e299,1.,'finite_physical_beta_covariance')
assert run('spatial_predict',[[0]],[1,2,1,0,0,1e100,0,1,0,0,0,1,1,0,0,1.2e-104,0,1.2e-104,5e-209])is None
proc.stdin.close();assert proc.wait()==0
src=ROOT/'crates/analytics/src/spatial.rs';test=ROOT/'crates/analytics/tests/candidates/spatial.rs';script=Path(__file__)
record={'schema_version':1,'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'PASS','layer':'isolated exact production spatial.rs via line-protocol Rust wrapper; public execute and native SQL are separate gates','seed':20261025,'calls':calls,'numeric_fields':fields,'maximum_scaled_error':maxerr,'groups':groups,'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'candidate_tests_sha256':hashlib.sha256(test.read_bytes()).hexdigest(),'harness_python':script.read_text(),'wrapper_source':Path(a.wrapper+'.rs').read_text(),'limitations':['supplied absolute covariance only; no hyperparameter estimation uncertainty','Isotropic four kernels in1..3 coordinate dimensions; no spatial model search or anisotropic parameter fitting','No installed-JMP comparison; no covariance-parameter SE or small-sample corrections']}
a.output.write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n');print(json.dumps({k:record[k]for k in ['status','calls','numeric_fields','maximum_scaled_error','groups']},indent=2))
