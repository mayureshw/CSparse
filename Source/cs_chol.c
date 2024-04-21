#include "cs.h"
#include "cholif.h"

static s_dword buf[CHOLBUFSZ];
static s_dword *p_buf0 = &buf[0], *p_buf1 = &buf[1], *p_buf;
static uint64_t cnt;

#define BUFCHK(CONTEXT) if ( cnt > CHOLBUFSZ ) { \
    printf("Buffer overfolw in %s\n",#CONTEXT); \
    exit(1); \
    }

/* L = chol (A, [pinv parent cp]), pinv is optional */
csn *cs_chol (const cs *A, const css *S)
{
    double *Cx ;
    csi top, i, p, k, n, *Li, *Lp, *cp, *pinv, *s, *c, *parent, *Cp, *Ci ;
    cs *L, *C, *E ;
    csn *N ;
    if (!CS_CSC (A) || !S || !S->cp || !S->parent) return (NULL) ;
    n = A->n ;
    N = cs_calloc (1, sizeof (csn)) ;       /* allocate result */
    c = cs_malloc (2*n, sizeof (csi)) ;     /* get csi workspace */
    cp = S->cp ; pinv = S->pinv ; parent = S->parent ;
    C = pinv ? cs_symperm (A, pinv, 1) : ((cs *) A) ;
    E = pinv ? C : NULL ;           /* E is alias for A, or a copy E=A(p,p) */
    if (!N || !c || !C) return (cs_ndone (N, E, c, 0, 0)) ;
    s = c + n ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    N->L = L = cs_spalloc (n, n, cp [n], 1, 0) ;    /* allocate result */
    if (!L) return (cs_ndone (N, E, c, 0, 0)) ;
    Lp = L->p ; Li = L->i ;
    for (k = 0 ; k < n ; k++) Lp [k] = c [k] = cp [k] ;

    p_buf0->l.u = n; // how many sends follow, also helps infer k
#ifdef CHOL_TRACE
    printf("n=%d\n",n);
    fflush(stdout);
#endif
    chol_send(p_buf0,1);

    for (k = 0 ; k < n ; k++)       /* compute L(k,:) for L*L' = C */
    {
        /* --- Nonzero pattern of L(k,:) ------------------------------------ */
        top = cs_ereach (C, k, parent, s, c) ;      /* find pattern of L(k,:) */
        
        p_buf = &buf[1]; // payload starts at index 1, will fill the header later
        cnt = 1;
        for (p = Cp [k] ; p < Cp [k+1] ; p++)       /* x = full(triu(C(:,k))) */
        {
            int cip = Ci[p];
            if ( cip <= k )
            {
                cnt++; BUFCHK(INITDX);
                p_buf->l.u = cip;
                p_buf->h.d = Cx[p];
                p_buf++;
#ifdef CHOL_TRACE
                printf("cip=%d\n",cip);
                printf("cxp=%.4e\n",Cx[p]);
                fflush(stdout);
#endif
            }
        }
        // even if 0, we have to send at least header to match count of n
        p_buf0->l.u = cnt - 1; // count excluding header
#ifdef CHOL_TRACE
        printf("initcnt=%d\n",cnt-1);
        fflush(stdout);
#endif

        /* --- Triangular solve --------------------------------------------- */
        csi c_k;
        c_k = c [k] ;
        cnt++; BUFCHK(TRSOLVE);
        p_buf->l.u = n - top;
        p_buf->h.u = c_k;
#ifdef CHOL_TRACE
        printf("n_minus_top=%d\n",n-top);
        printf("c_k=%d\n",c_k);
        fflush(stdout);
#endif
        chol_send(p_buf0,cnt);
        for ( ; top < n ; top++)    /* solve L(0:k-1,0:k-1) * x = C(:,k) */
        {
            p_buf = p_buf0; cnt = 0;
            i = s [top] ;               /* s [top..n-1] is pattern of L(k,:) */
            int c_i = c[i];

            cnt++; BUFCHK(TRSOLVE);
            p_buf->l.u = i;
            p_buf->h.u = c_i;
            p_buf++;
            cnt++; BUFCHK(TRSOLVE);
            p_buf->l.u = Lp[i];
            p_buf++;
#ifdef CHOL_TRACE
        printf("i=%d\n",i);
        printf("c_i=%d\n",c_i);
        printf("Lp_i=%d\n",Lp[i]);
        fflush(stdout);
#endif

            Li [c_i] = k ;                /* store L(k,i) in column i */
            char use_l = 1; // a toggle flag to fill l or h
            for (p = Lp [i] + 1 ; p < c_i ; p++)
            {
                if ( use_l )
                {
                    cnt++; BUFCHK(TRSOLVE);
                    p_buf->l.u = Li[p];
                    use_l = 0;
                }
                else
                {
                    p_buf->h.u = Li[p];
                    p_buf++;
                    use_l = 1;
                }
#ifdef CHOL_TRACE
                printf("Li_p=%d\n",Li[p]);
                fflush(stdout);
#endif
            }
            chol_send(p_buf0,cnt);
            c [i]++ ;
        }
        /* --- Compute L(k,k) ----------------------------------------------- */
        Li [c_k] = k ;                /* store L(k,k) = sqrt (d) in column k */
        c[k]++;
    }
    chol_flush();
    chol_read(L->x,cp[n]);
    Lp [n] = cp [n] ;               /* finalize L */
    return (cs_ndone (N, E, c, 0, 1)) ; /* success: free E,s,x; return N */
}
