"""Seeded Python reference experiment for the archived R equations, not an R rerun."""
from pathlib import Path
import json,argparse
import numpy as np

def network(beta,h,rng,nc=15):
 n=len(beta);p=nc*np.outer(beta,beta)/(beta.mean()**2*(n-1));p=np.clip(p,0,1);p[h[:,None]==h[None,:]]=0
 A=np.triu(rng.random((n,n))<p,1);A=A|A.T
 assert np.array_equal(A,A.T) and not np.any(np.diag(A)) and not np.any(A[h[:,None]==h[None,:]])
 return A

def simulate(beta,h,A,rng,alpha=(.1,.01,.01),nt=100):
 n=len(beta);x=np.zeros(n,dtype=np.int8);x[rng.choice(n,round(.005*n),replace=False)]=2;counts=[np.bincount(x,minlength=4)]
 for _ in range(1,nt):
  new=x.copy();S=np.flatnonzero(x==0);E=np.flatnonzero(x==1);I=np.flatnonzero(x==2)
  new[E[rng.random(len(E))<.4]]=2;new[I[rng.random(len(I))<.2]]=3
  household_counts=np.bincount(h[I],minlength=int(h.max())+1)
  ph=1-(1-alpha[0])**household_counts[h[S]];pc=1-(1-alpha[1])**A[np.ix_(S,I)].sum(axis=1)
  M=alpha[2]*15*np.outer(beta[S],beta[I])/(beta.mean()**2*(n-1))
  assert np.all((M>=0)&(M<1))
  pr=-np.expm1(np.log1p(-M).sum(axis=1));p=1-(1-ph)*(1-pc)*(1-pr)
  assert np.all((p>=0)&(p<=1));new[S[rng.random(len(S))<p]]=1;x=new;counts.append(np.bincount(x,minlength=4))
 counts=np.asarray(counts);assert np.all(counts.sum(axis=1)==n) and np.all(np.diff(counts[:,0])<=0) and np.all(np.diff(counts[:,3])>=0)
 return {'peak_infectious':int(counts[:,2].max()),'peak_day':int(counts[:,2].argmax()),'ever_infected_by_day99':int(n-counts[-1,0]),'infectious_day99':int(counts[-1,2])}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--replicates',type=int,default=50);args=parser.parse_args();cases={str(i):[] for i in range(1,5)}
 for seed in range(args.replicates):
  rng=np.random.default_rng(20261101+seed);h=np.repeat(np.arange(1000),rng.integers(1,6,1000))[:1000];b=rng.uniform(0,1,1000)
  for i,(beta,alpha) in enumerate([(b,(.1,.01,.01)),(b,(0,0,.04)),(np.repeat(b.mean(),1000),(.1,.01,.01)),(np.repeat(b.mean(),1000),(0,0,.04))],start=1):
   A=network(beta,h,rng);cases[str(i)].append(simulate(beta,h,A,rng,alpha))
 summary={}
 for case,rows in cases.items():
  summary[case]={}
  for key in rows[0]:
   vals=np.array([x[key] for x in rows]);m=float(vals.mean());d=float(1.96*vals.std(ddof=1)/np.sqrt(len(vals)));summary[case][key]={'mean':m,'ci95_mean':[m-d,m+d]}
 result={'reviewed_on':'2026-11-01','scope':'New seeded Python reference simulations of the archived equations; original R runs not reproduced.','population':1000,'days':100,'initial_infectious':5,'replicates_per_case':args.replicates,'master_seed':20261101,'cases':summary,'invariants_passed_every_run':['population conserved','susceptibles nonincreasing','removed nondecreasing','network symmetric and no self/household links','transmission probabilities valid'],'interpretation':'Conditional simulation outcomes only. Removing household/network transmission also raises alpha_random from 0.01 to 0.04, so this is not a one-factor causal test of contact structure. Replication intervals concern Monte Carlo mean precision, not real-world epidemic uncertainty.'}
 out=Path(__file__).resolve().parent/'REVIEW_VERIFICATION.json';out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
