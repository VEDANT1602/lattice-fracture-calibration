"""Single-threaded repeated solves and certified scalar queries on a cut square network."""
import os
for variable in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[variable]="1"
import argparse
import csv
import json
import platform
import sys
from pathlib import Path
from time import perf_counter

import numpy as np
import scipy
from scipy.ndimage import gaussian_filter
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu

ROOT=Path(__file__).resolve().parent
FAMILIES=("fixed_local","distributed_mild","distributed_strong")


def write_csv(path,rows):
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


class Network:
    def __init__(self,n):
        if n % 2: raise ValueError("use an even width")
        self.n=n
        self.dim=n*(n-2)
        idx=np.arange(n*n).reshape(n,n)
        a=np.r_[idx[:,:-1].ravel(),idx[:-1,:].ravel()]
        b=np.r_[idx[:,1:].ravel(),idx[1:,:].ravel()]
        x=(a%n+b%n)/2
        y=(a//n+b//n)/2
        at_cut=(y==(n-1)/2)&(x==n//4)
        self.h=np.zeros(self.dim)
        self.h[[int(a[at_cut].item())-n,int(b[at_cut].item())-n]]=[-1,1]
        keep=~((y==(n-1)/2)&(x<=n//4))
        self.a,self.b=a[keep]-n,b[keep]-n
        self.x,self.y=x[keep],y[keep]
        self.m=len(self.a)
        near=np.hypot(self.x-n//4,self.y-(n-1)/2)
        self.local=np.lexsort((np.arange(self.m),near))[:16]
        rows,cols,refs,signs=[],[],[],[]
        for i,j,s in [(self.a,self.a,1),(self.b,self.b,1),(self.a,self.b,-1),(self.b,self.a,-1)]:
            valid=(i>=0)&(j>=0)&(i<self.dim)&(j<self.dim)
            rows.append(i[valid]);cols.append(j[valid]);refs.append(np.flatnonzero(valid))
            signs.append(np.full(valid.sum(),s))
        keys,inv=np.unique(np.concatenate(cols)*self.dim+np.concatenate(rows),return_inverse=True)
        self.ptr=np.r_[0,np.cumsum(np.bincount(keys//self.dim,minlength=self.dim))].astype(np.int32)
        self.indices=(keys%self.dim).astype(np.int32)
        self.refs=np.concatenate(refs);self.signs=np.concatenate(signs);self.inv=inv
        self.ones=np.ones(self.m)

    def matrix(self,k):
        values=np.bincount(self.inv,weights=k[self.refs]*self.signs)
        return csc_matrix((values,self.indices,self.ptr),shape=(self.dim,self.dim))

    def factor_base(self):
        self.base=splu(self.matrix(self.ones))
        self.x0=self.base.solve(self.h)
        self.r0=float(self.h@self.x0)

    def gradients(self,x):
        full=np.pad(x,(self.n,self.n))
        return full[self.b+self.n]-full[self.a+self.n]

    def setup_screen(self):
        v=self.gradients(self.x0)
        self.w=v*v/self.r0
        np.testing.assert_allclose(self.w.sum(),1,atol=2e-12)

    def setup_woodbury(self):
        # A single fixed support is known from the stochastic input description.
        i,j=self.a[self.local],self.b[self.local]
        assert (i>=0).all() and (j<self.dim).all()
        columns=np.zeros((self.dim,len(i)),order="F")
        columns[i,np.arange(len(i))]=-1
        columns[j,np.arange(len(j))]=1
        inverse_times_B=self.base.solve(columns)
        self.schur=inverse_times_B[j,:]-inverse_times_B[i,:]
        self.f=self.x0[j]-self.x0[i]

    def field(self,family,seed):
        rng=np.random.default_rng(seed)
        if family=="fixed_local":
            k=self.ones.copy()
            k[self.local]=np.exp(rng.uniform(-np.log(4),np.log(4),16))
        else:
            lattice=gaussian_filter(rng.standard_normal((self.n,self.n)),sigma=3.0,mode="reflect").ravel()
            edge=np.sign(lattice[self.a+self.n]+lattice[self.b+self.n])
            ratio=1.25 if family=="distributed_mild" else 4.0
            k=np.power(ratio,edge)
        return k

    def screen(self,k,local):
        if local:
            v=k[self.local];w=self.w[self.local]
            a=1+w@(v-1);b=1+w@(1/v-1)
        else:
            a=self.w@k;b=self.w@(1/k)
        return (1+self.r0/a)/2,(1+self.r0*b)/2

    def exact(self,k):
        return (1+self.h@splu(self.matrix(k)).solve(self.h))/2

    def woodbury(self,k):
        delta=k[self.local]-1
        small=np.eye(len(delta))+delta[:,None]*self.schur
        change=self.f@np.linalg.solve(small,delta*self.f)
        return (1+self.r0-change)/2

    def pcg(self,k,mode):
        A=self.matrix(k).tocsr()
        x=self.x0.copy()
        kmin=min(1,float(k.min()))
        r=self.h-A@x
        z=self.base.solve(r)
        p=z.copy()
        rho=float(r@z)
        hnorm=np.linalg.norm(self.h)
        for step in range(151):
            variational=float(self.h@x+r@x)
            lo=(1+variational)/2
            hi=lo+rho/(2*kmin)
            accept=mode=="output" and lo>0 and hi<=lo*1.02**2
            if mode=="tail":
                threshold=(1+self.r0)/2
                accept=(hi<threshold-1e-11 or lo>threshold+1e-11)
            if accept or np.linalg.norm(r)<=1e-10*hnorm:
                # Refresh a true residual before recording the certificate.
                r=self.h-A@x
                rr=float(r@self.base.solve(r))
                lo=(1+self.h@x+r@x)/2
                hi=lo+rr/(2*kmin)
                return np.sqrt(lo*hi),lo,hi,step
            q=A@p
            alpha=rho/float(p@q)
            x+=alpha*p;r-=alpha*q
            z=self.base.solve(r)
            newrho=float(r@z)
            p=z+(newrho/rho)*p;rho=newrho
        raise RuntimeError("PCG did not converge")


def benchmark(n,nsamples,out):
    t=perf_counter();g=Network(n);geom_ms=(perf_counter()-t)*1000
    t=perf_counter();g.factor_base();base_ms=(perf_counter()-t)*1000
    t=perf_counter();g.setup_screen();screen_ms=(perf_counter()-t)*1000
    t=perf_counter();g.setup_woodbury();woodbury_ms=(perf_counter()-t)*1000
    ref=(1+g.r0)/2
    setup=dict(n=n,unknowns=g.dim,edges=g.m,geometry_ms=geom_ms,uniform_factor_solve_ms=base_ms,
               screen_extra_ms=screen_ms,woodbury_extra_ms=woodbury_ms,uniform_output=ref)
    rows=[]
    tasks=[(family,i) for family in FAMILIES for i in range(nsamples)]
    rng_order=np.random.default_rng(93313)
    for case in rng_order.permutation(len(tasks)):
        family,i=tasks[case];seed=81000+i
        k=g.field(family,seed);local=(family=="fixed_local")
        row=dict(n=n,family=family,seed=seed)
        methods=['direct','screen','tight_pcg','cert_pcg','tail_pcg']
        if local:methods.append('woodbury')
        for method in rng_order.permutation(methods):
            t=perf_counter()
            if method=='direct':
                val=g.exact(k)
            elif method=='screen':
                lo,hi=g.screen(k,local);val=np.sqrt(lo*hi)
                row['screen_lower']=lo;row['screen_upper']=hi
            elif method=='woodbury':
                val=g.woodbury(k)
            else:
                mode={'cert_pcg':'output','tight_pcg':'tight','tail_pcg':'tail'}[method]
                val,lo,hi,it=g.pcg(k,mode)
                row[method+'_lower']=lo;row[method+'_upper']=hi;row[method+'_iterations']=it
            elapsed=1000*(perf_counter()-t)
            row[method]=val;row[method+'_ms']=elapsed
        if not local:row['woodbury']=row['woodbury_ms']=''
        row['screen_pass']=int(row['screen_upper']<=row['screen_lower']*1.02**2)
        row['screen_tail_pass']=int(row['screen_lower']>ref+1e-11 or row['screen_upper']<ref-1e-11)
        for method,stop in [('screen','screen_pass'),('screen_tail','screen_tail_pass')]:
            pcg='cert_pcg' if method=='screen' else 'tail_pcg'
            row[method+'_hybrid_ms']=row['screen_ms']+(0 if row[stop] else row[pcg+'_ms'])
        row['hybrid_value']=row['screen'] if row['screen_pass'] else row['cert_pcg']
        if not(row['screen_lower']-1e-9 <= row['direct'] <= row['screen_upper']+1e-9):
            raise AssertionError('screen bracket')
        for p in ['tight_pcg','cert_pcg','tail_pcg']:
            if not row[p+'_lower']-1e-9<=row['direct']<=row[p+'_upper']+1e-9:
                raise AssertionError('PCG bracket')
        assert abs(row['tight_pcg']/row['direct']-1)<1e-8
        assert abs(row['hybrid_value']/row['direct']-1)<.020000001
        if local: assert abs(row['woodbury']/row['direct']-1)<1e-10
        if abs(row['direct']-ref)>1e-10: assert (row['tail_pcg']>ref)==(row['direct']>ref)
        rows.append(row)
    rows.sort(key=lambda r:(r['family'],r['seed']))
    columns=sorted(rows[0])
    write_csv(out/f'queries_{n}.csv',[{c:r[c] for c in columns} for r in rows])
    return setup


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--sizes',nargs='+',type=int,default=[64,208])
    cli.add_argument('--samples',type=int,default=48)
    cli.add_argument('--out',type=Path,default=ROOT)
    opts=cli.parse_args()
    setups=[]
    for n in opts.sizes:
        setups.append(benchmark(n,opts.samples,opts.out))
        write_csv(opts.out/'setup.csv',setups)
        print('Completed',n,flush=True)
    (opts.out/'environment.json').write_text(json.dumps(dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),threads=1),indent=2)+'\n')
