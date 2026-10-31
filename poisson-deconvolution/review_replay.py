"""Independent Python reconstruction of the archived Poisson objective; not an R rerun."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np,pandas as pd
from scipy.interpolate import BSpline
from scipy.stats import lognorm
from scipy.optimize import minimize
from scipy.linalg import solve

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
 a=pd.read_csv(args.input,sep=r'\s+');y=a.nhs.to_numpy(float);t=a.julian.to_numpy(int);n=len(y);K=80
 assert np.all(np.diff(t)==1) and np.all(y>=0)
 time=np.arange(t.min()-30,t.max()+1);span=np.ptp(time);knots=np.linspace(time.min()-.01*span,time.max()+.01*span,K+4)
 B=np.column_stack([np.nan_to_num(BSpline.basis_element(knots[k:k+5],extrapolate=False)(time),nan=0) for k in range(K)])
 pdur=lognorm.pdf(np.arange(1,81),s=.469,scale=np.exp(3.151));pdur/=pdur.sum()
 X=np.array([np.sum(B[ti-np.arange(1,min(30+i,80)+1)-time[0]]*pdur[:min(30+i,80),None],axis=0) for i,ti in enumerate(t)])
 D=np.diff(np.eye(K),n=2,axis=0);S=D.T@D
 def value_grad(g,lam,w):
  b=np.exp(g);m=np.maximum(X@b,1e-9);pen=S@b
  return float(np.sum(w*(m-y*np.log(m)))+.5*lam*b@pen),b*(X.T@(w*(1-y/m))+lam*pen)
 def fit(lam,start,w):return minimize(lambda g:value_grad(g,lam,w),start,jac=True,method='BFGS',options={'maxiter':2000,'gtol':1e-5})
 w=np.ones(n);g=np.repeat(np.log(y.mean()),K);eps=1e-5
 analytic=value_grad(g,5e-5,w)[1];fd=np.array([(value_grad(g+np.eye(K)[j]*eps,5e-5,w)[0]-value_grad(g-np.eye(K)[j]*eps,5e-5,w)[0])/(2*eps) for j in range(K)])
 initial=fit(5e-5,g,w);grid=[];solutions=[]
 for loglam in np.linspace(-13,-7,50):
  lam=np.exp(loglam);f=fit(lam,initial.x,w);b=np.exp(f.x);m=X@b;H=X.T@((y/m**2)[:,None]*X);edf=float(np.trace(solve(H+lam*S,H,assume_a='pos')));nll=float(np.sum(m-y*np.log(m)));bic=2*nll+np.log(n)*edf
  grid.append({'log_lambda':float(loglam),'lambda':float(lam),'BIC_without_poisson_constants':float(bic),'edf':edf,'optimizer_success':bool(f.success),'gradient_inf':float(np.max(abs(f.jac))),'iterations':int(f.nit)});solutions.append(f)
 # Numerical convergence is assessed by gradient as well as SciPy's termination flag.
 eligible=[i for i,z in enumerate(grid) if z['gradient_inf']<1e-3];assert eligible,'No adequately stationary candidate'
 best=min(eligible,key=lambda i:grid[i]['BIC_without_poisson_constants']);f=solutions[best];lam=grid[best]['lambda'];curve=B@np.exp(f.x);mu=X@np.exp(f.x)
 rng=np.random.default_rng(20261101);boots=[];boot_grad=[]
 for _ in range(200):
  weights=np.bincount(rng.integers(0,n,n),minlength=n);bf=fit(lam,f.x,weights);boot_grad.append(float(np.max(abs(bf.jac))));boots.append(B@np.exp(bf.x))
 boots=np.asarray(boots);ci=np.quantile(boots,[.025,.975],axis=0)
 args.out.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame({'day':time,'fatal_infections':curve,'lower':ci[0],'upper':ci[1]}).to_csv(args.out.parent/'poisson_curve_private.csv',index=False)
 pd.DataFrame(grid).to_csv(args.out.parent/'poisson_selection.csv',index=False)
 result={'reviewed_on':'2026-11-01','scope':'Independent Python implementation of the archived R model on the recovered input; original R script not executed in this check.','input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),'observations':n,'first_day':int(t.min()),'last_day':int(t.max()),'deaths_total':int(y.sum()),'basis_coefficients':K,'infection_grid_points':len(time),'delay_probability_sum':float(pdur.sum()),'central_difference_max_gradient_error':float(np.max(abs(analytic-fd))),'selected':grid[best],'eligible_gradient_below_1e_minus3_count':len(eligible),'candidate_count':len(grid),'fatal_infection_peak_day':int(time[np.argmax(curve)]),'fitted_death_peak_day':int(t[np.argmax(mu)]),'peak_lag_days':int(t[np.argmax(mu)]-time[np.argmax(curve)]),'in_sample_death_RMSE':float(np.sqrt(np.mean((mu-y)**2))),'bootstrap_replicates':200,'bootstrap_seed':20261101,'bootstrap_max_gradient':float(max(boot_grad)),'bootstrap_fits_gradient_above_1e_minus3':int(sum(x>=1e-3 for x in boot_grad)),'uncertainty':'Fixed-smoothing weighted-case bootstrap; no simultaneous, selection, delay-kernel or dependence uncertainty.','no_predictive_holdout':True}
 args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
