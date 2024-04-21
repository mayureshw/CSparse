#include "cs.h"
#include "cholif.h"

#define WRBUFSZ (1<<26)

static s_dword *buf, *p_first, *p_last, *p_buf;

// To be invoked after writing 1 record to *p_buf
void bufnext()
{
    if ( p_buf == p_last )
    {
        chol_send( p_first, WRBUFSZ );
        p_buf = p_first;
    }
    else p_buf++;
}

// Expects the pointer to be at a position that is not yet written to
// This helps decide whether flush is needed at all
// For this, ensure to do bufnext after populating each dword
void bufflush()
{
    if ( p_buf != p_first ) chol_send( p_first, p_buf - p_first );
}

void initbuf()
{
    buf = malloc(WRBUFSZ*sizeof(s_dword));
    if ( !buf )
    {
        printf("Could not initialize buf for size %d\n",WRBUFSZ);
        exit(1);
    }
    p_buf = p_first = &buf[0];
    p_last = &buf[WRBUFSZ-1];
}

/* L = chol (A, [pinv parent cp]), pinv is optional */
csn *cs_chol (const cs *A, const css *S)
{
    initbuf();
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

    p_buf->l.u = n; // how many sends follow, also helps infer k
#ifdef CHOL_TRACE
    printf("n=%d\n",n);
    fflush(stdout);
#endif
    bufnext();

    for (k = 0 ; k < n ; k++)       /* compute L(k,:) for L*L' = C */
    {
        /* --- Nonzero pattern of L(k,:) ------------------------------------ */
        top = cs_ereach (C, k, parent, s, c) ;      /* find pattern of L(k,:) */
        
        // TODO: With bufnext design, can't wait to populate header later
        // So have to walk over Cp twice to get the count first. Can be improved.
        uint64_t initcnt = 0;
        for (p = Cp [k] ; p < Cp [k+1] ; p++) if ( Ci[p] <=k ) initcnt++;
#ifdef CHOL_TRACE
        printf("initcnt=%d\n",initcnt);
        fflush(stdout);
#endif
        p_buf->l.u = initcnt;
        bufnext();

        for (p = Cp [k] ; p < Cp [k+1] ; p++)       /* x = full(triu(C(:,k))) */
        {
            int cip = Ci[p];
            if ( cip <= k )
            {
                p_buf->l.u = cip;
                p_buf->h.d = Cx[p];
                bufnext();
#ifdef CHOL_TRACE
                printf("cip=%d\n",cip);
                printf("cxp=%.4e\n",Cx[p]);
                fflush(stdout);
#endif
            }
        }

        /* --- Triangular solve --------------------------------------------- */
        csi c_k;
        c_k = c [k] ;
        p_buf->l.u = n - top;
        p_buf->h.u = c_k;
        bufnext();
#ifdef CHOL_TRACE
        printf("n_minus_top=%d\n",n-top);
        printf("c_k=%d\n",c_k);
        fflush(stdout);
#endif
        for ( ; top < n ; top++)    /* solve L(0:k-1,0:k-1) * x = C(:,k) */
        {
            i = s [top] ;               /* s [top..n-1] is pattern of L(k,:) */
            int c_i = c[i];

            p_buf->l.u = i;
            p_buf->h.u = c_i;
            bufnext();

            p_buf->l.u = Lp[i];
            bufnext();
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
                    p_buf->l.u = Li[p];
                    use_l = 0;
                }
                else
                {
                    p_buf->h.u = Li[p];
                    bufnext();
                    use_l = 1;
                }
#ifdef CHOL_TRACE
                printf("Li_p=%d\n",Li[p]);
                fflush(stdout);
#endif
            }
            if ( ! use_l ) bufnext();
            c [i]++ ;
        }
        /* --- Compute L(k,k) ----------------------------------------------- */
        Li [c_k] = k ;                /* store L(k,k) = sqrt (d) in column k */
        c[k]++;
    }
    bufflush();
    chol_read(L->x,cp[n]);
    Lp [n] = cp [n] ;               /* finalize L */
    return (cs_ndone (N, E, c, 0, 1)) ; /* success: free E,s,x; return N */
}
