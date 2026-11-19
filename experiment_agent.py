import itertools,random,csv,json,time,math,sys
from collections import defaultdict

def rref(A,p):
 A=[row[:] for row in A]; m=len(A); n=len(A[0]) if m else 0; piv=[]; ops=0
 for j in range(n):
  k=next((k for k in range(len(piv),m) if A[k][j]%p),None)
  if k is None: continue
  i=len(piv); A[i],A[k]=A[k],A[i]; z=pow(A[i][j],-1,p); ops+=1
  for c in range(j,n): A[i][c]=A[i][c]*z%p; ops+=1
  for k in range(m):
   if k!=i and A[k][j]:
    t=A[k][j]
    for c in range(j,n): A[k][c]=(A[k][c]-t*A[i][c])%p; ops+=2
  piv.append(j)
  if len(piv)==m: break
 return A,piv,ops

def kernel(A,n,p):
 if not A:return [[int(i==j) for i in range(n)] for j in range(n)],0
 R,piv,ops=rref(A,p); free=[j for j in range(n) if j not in piv]; B=[]
 for f in free:
  v=[0]*n;v[f]=1
  for i,j in enumerate(piv):v[j]=-R[i][f]%p
  B.append(v)
 return B,ops

def monomials(d,b):
 if d==1: yield (b,);return
 for v in range(b+1):
  for w in monomials(d-1,b-v):yield (v,)+w

def match(T,r,mask):
 for perm in itertools.permutations(range(r)):
  if all((j==i if j<r else (i,j-r) not in mask) for i,j in zip(perm,T)):return True
 return False

def prod(A,B,p):return [[sum(x*y for x,y in zip(row,col))%p for col in zip(*B)] for row in A]

def trial(p,m,n,K,r,h,seed,maxb=3):
 rng=random.Random(seed//100); mask_rng=random.Random(seed); mask=set(mask_rng.sample(list(itertools.product(range(r),range(n-r))),h));
 C=[[int(i==j) if j<r else rng.randrange(p) for j in range(n)] for i in range(r)]
 S=[[rng.randrange(p) for _ in range(r)] for _ in range(m)]; E=prod(S,C,p)
 # force full rank S; avoid trivial x planted in a basis by selecting random x* with x*_0=1
 while len(rref(S,p)[1])<r:S=[[rng.randrange(p) for _ in range(r)] for _ in range(m)]
 E=prod(S,C,p)
 xstar=[1]+[rng.randrange(p) for _ in range(K-1)]; mats=[None]+[[[rng.randrange(p) for _ in range(n)] for _ in range(m)] for _ in range(K-1)]
 mats[0]=[[(E[i][j]-sum(xstar[k]*mats[k][i][j] for k in range(1,K)))%p for j in range(n)] for i in range(m)]
 return solve(p,mats,r,{(i,j):C[i][r+j] for i,j in mask},maxb,seed)

def solve(p,mats,r,hints,maxb=None,seed=None):
 """Solve supplied public matrices with hints keyed by zero-based (row, nonidentity-column)."""
 start=time.perf_counter()
 if p<2 or any(p%i==0 for i in range(2,math.isqrt(p)+1)):raise ValueError("q must be prime")
 if maxb is None:maxb=min(3,p-1)
 if not 1<=maxb<=p-1:raise ValueError("max_b must satisfy 1 <= max_b < q")
 K=len(mats);m=len(mats[0]);n=len(mats[0][0]);mask=set(hints);h=len(mask)
 if not 1<=r<=min(m,n) or any(len(M)!=m or any(len(row)!=n for row in M) for M in mats):raise ValueError("matrix dimensions or r invalid")
 if any(not (0<=i<r and 0<=j<n-r and 0<=v<p) for (i,j),v in hints.items()):raise ValueError("invalid hint")
 C=[[int(i==j) if j<r else hints.get((i,j-r),0) for j in range(n)] for i in range(r)]
 full=[j for j in range(n-r) if all((i,j) in mask for i in range(r))]; rem=[j for j in range(n-r) if j not in full]; keep=list(range(r))+[r+j for j in rem]
 linear_start=time.perf_counter()
 linear=[[ (M[i][r+j]-sum(C[t][r+j]*M[i][t] for t in range(r)))%p for M in mats] for j in full for i in range(m)]
 basis,op1=kernel(linear,K,p);d=len(basis);linear_seconds=time.perf_counter()-linear_start
 if not d:return {'status':'no_nonzero_solution_after_linear_constraints','params':[p,m,n,K,r,h],'linear_ops':op1,'linear_seconds':linear_seconds}
 if len(keep)==r:
  # All nonidentity columns are known. Search the reduced projective space.
  # A valid candidate already satisfies the known-column relations.
  checked=0
  for lead in range(d):
   for tail in itertools.product(range(p),repeat=d-lead-1):
    u=[0]*lead+[1]+list(tail);X=[sum(u[t]*basis[t][j] for t in range(d))%p for j in range(K)]
    Em=[[sum(X[k]*mats[k][i][j] for k in range(K))%p for j in range(n)] for i in range(m)]
    checked+=1
    if len(rref(Em,p)[1])==r:return {'status':'solved_by_full_column_reduction','params':[p,m,n,K,r,h],'d':d,'linear_ops':op1,'projective_candidates_checked':checked,'solution':{'x':X,'valid':True}}
    if checked>=10000:return {'status':'candidate_search_cap','params':[p,m,n,K,r,h],'d':d,'linear_ops':op1,'projective_candidates_checked':checked}
  return {'status':'no_rank_r_solution_in_chart','params':[p,m,n,K,r,h],'d':d,'linear_ops':op1,'projective_candidates_checked':checked}
 # transformed matrices and normalization, count straightforward multiplication and subtraction
 transform_start=time.perf_counter(); reduced=[];ops2=0
 for v in basis:
  mat=[[sum(v[k]*mats[k][i][j] for k in range(K))%p for j in keep] for i in range(m)];ops2+=m*len(keep)*2*K
  for i in range(m):
   for z,j in enumerate(rem):
    for t in range(r):
     if (t,j) in mask:mat[i][r+z]=(mat[i][r+z]-C[t][r+j]*mat[i][t])%p;ops2+=2
  reduced.append(mat)
 transform_seconds=time.perf_counter()-transform_start
 matching_start=time.perf_counter()
 n0=len(keep);mask0={(i,rem.index(j)) for i,j in mask if j in rem};Ts=[T for T in itertools.combinations(range(n0),r) if match(T,r,mask0)];Tmap={T:i for i,T in enumerate(Ts)};Ls=[L for L in itertools.combinations(range(n0),r+1) if any(L[:i]+L[i+1:] in Tmap for i in range(r+1))];eq=[]
 matching_seconds=time.perf_counter()-matching_start
 equation_start=time.perf_counter()
 for row in range(m):
  for L in Ls:
   e={}
   for i,j in enumerate(L):
    T=L[:i]+L[i+1:]
    if T not in Tmap:continue
    for k in range(d):
     a=((-1)**i*reduced[k][row][j])%p
     if a:e[k,Tmap[T]]=a
   if e:eq.append(e)
 equation_seconds=time.perf_counter()-equation_start
 records=[]
 for b in range(1,min(maxb,p-1)+1):
  build_start=time.perf_counter()
  mon=list(monomials(d,b)); lower=list(monomials(d,b-1));mid={v:i for i,v in enumerate(mon)};N=len(mon)*len(Ts);A=[]
  if N*len(lower)*len(eq)>2000000:
   records.append({"b":b,"status":"resource_cap", "rows":len(lower)*len(eq),"cols":N});break
  for v in lower:
   for e in eq:
    row=[0]*N
    for (k,t),coeff in e.items():
     w=list(v);w[k]+=1;row[mid[tuple(w)]*len(Ts)+t]=coeff
    A.append(row)
  build_seconds=time.perf_counter()-build_start;solve_start=time.perf_counter()
  ker,ops=kernel(A,N,p);solve_seconds=time.perf_counter()-solve_start;rank=N-len(ker);rec={'b':b,'rows':len(A),'cols':N,'rank':rank,'nullity':len(ker),'elim_ops':ops,'inversions':rank,'build_seconds':build_seconds,'solve_seconds':solve_seconds,'nnz':sum(x!=0 for row in A for x in row)};records.append(rec)
  if len(ker)==1:
   v=ker[0]; yindex=Tmap[tuple(range(r))];sol=None
   for k in range(d):
    wk=[0]*d;wk[k]=b;z=v[mid[tuple(wk)]*len(Ts)+yindex]
    if z:
     xs=[]
     for j in range(d):
      w=[0]*d;w[k]=b-1;w[j]+=1;xs.append(v[mid[tuple(w)]*len(Ts)+yindex]*pow(z,-1,p)%p)
     X=[sum(xs[t]*basis[t][j] for t in range(d))%p for j in range(K)]
     Em=[[sum(X[k]*mats[k][i][j] for k in range(K))%p for j in range(n)] for i in range(m)]
         # check hints in normalized chart through S=first r columns, solve S*C entries using planted hint and unknown values
     # when only partially hinted, find the first independent rows of S and recover C
     St=[[Em[i][j] for j in range(r)] for i in range(m)]; chosen=next((inds for inds in itertools.combinations(range(m),r) if len(rref([St[i] for i in inds],p)[1])==r),None)
     valid=False
     if chosen:
      invR,piv,_=rref([[St[i][j] for j in range(r)]+[int(a==b) for b in range(r)] for a,i in enumerate(chosen)],p)
      inv=[row[r:] for row in invR[:r]]; recovered=prod(inv,[Em[i] for i in chosen],p); valid=len(rref(Em,p)[1])==r and prod(St,recovered,p)==Em and all(recovered[i][r+j]==C[i][r+j] for i,j in mask)
     sol={'x':X,'valid':valid};break
   rec['solution']=sol
   if sol and sol['valid']:break
 return {'params':[p,m,n,K,r,h],'seed':seed,'hint_values':[[i,j,hints[i,j]] for i,j in sorted(hints)],'seconds':{'linear':linear_seconds,'transform':transform_seconds,'matching':matching_seconds,'equations':equation_seconds,'total':time.perf_counter()-start},'mask':sorted([list(x) for x in mask]),'full':full,'d':d,'n0':n0,'s':len(Ts),'blocks':len(Ls),'linear_ops':op1,'transform_ops':ops2,'matching_tests':math.comb(n0,r),'macaulay':records}

if __name__=='__main__':
 if len(sys.argv)>1:
  with open(sys.argv[1]) as f:data=json.load(f)
  hints={(int(a),int(b)):int(v) for a,b,v in data['hints']}
  if len(hints)!=len(data['hints']):raise ValueError('duplicate hint positions')
  print(json.dumps(solve(data['q'],data['matrices'],data['r'],hints,data.get('max_b',min(3,data['q']-1))),indent=2));sys.exit(0)
 # Candidate sets; all arithmetic in prime fields, rows m, columns n, K public matrices, target r
 configs=[(7,4,5,5,2,2),(7,4,5,6,2,3),(7,5,6,7,2,3),(7,5,6,8,2,4),(7,5,6,7,3,5)]
 results=[]
 for z,cfg in enumerate(configs):
  for t in range(10):
   t0=time.time(); a=trial(*cfg,seed=10000*z+t);results.append(a);print(z,t,a['d'],[(q['b'],q['cols'],q['nullity']) for q in a['macaulay']],round(time.time()-t0,2),flush=True)
 with open('experiment_results_agent.json','w') as f:json.dump(results,f,indent=2)
