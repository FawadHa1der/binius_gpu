#include "prod_check.h"
#include <math.h>   
#include "../debug_utils.h"
#include "univariate_poly.h"
#include "compressed_poly.h"
// A "default" constructor => sets N=0, arrays= null, claim= zero, ...
ProdCheck prodcheck_default(void)
{
    ProdCheck pc;
    pc.N=0;
    pc.p_polys=NULL;
    pc.q_polys=NULL;
    pc.claim= f128_zero();
    // points_init(&pc.challenges);
    pc.num_vars=0;
    pc.has_cached=0; 
    return pc;
}


// Helper: compute integer log2
static inline size_t int_log2(size_t x){
    // For simplicity, use standard log2 and cast
    // must be exact for 2^n
    double lf= log2((double)x);
    size_t r= (size_t)(lf+0.5);
    return r;
}


ProdCheck prodcheck_new(MLE_POLY* p_arr, 
    MLE_POLY* q_arr,
                        size_t N,
                        F128 claim,
                        int check_init_claim)
{
    // Check they all have same length for the first pair
    size_t len0= p_arr[0].len;
    size_t num_vars= int_log2(len0);

    // confirm each p_i, q_i has length= 2^num_vars
    for(size_t i=0;i<N;i++){
        if(p_arr[i].len != (1ULL<<num_vars) || q_arr[i].len!=(1ULL<<num_vars)){
            fprintf(stderr,"All polynomials must have size 2^num_vars\n");
            exit(1);
        }
    }

    // Optionally verify claim
    if(check_init_claim){
        // compute expected by sum( p[i][j]* q[i][j] )
        F128 expected= f128_zero();
        for(size_t i=0;i<N;i++){
            F128 partial= f128_zero();
            for(size_t j=0;j< p_arr[i].len;j++){

                F128 pm= p_arr[i].coeffs[j];//  mlp_get(&(p_arr[i]), j);
                F128 qm= q_arr[i].coeffs[j]; // mlp_get(&(q_arr[i]), j);
                partial= f128_add(partial, f128_mul(pm, qm));
            }
            expected= f128_add(expected, partial);
        }
        if(!f128_eq(expected, claim)){
            // print them
            f128_print("expected", expected);
            f128_print( "claim", claim);
            fprintf(stderr,"Initial claim does not match computed value\n");
            exit(1);
        }
    }

    ProdCheck pc= prodcheck_default();
    pc.N= N;
    pc.challenges= points_init(num_vars); 
    // allocate new arrays, copy references
    pc.p_polys= (MLE_POLY*)malloc(sizeof(MLE_POLY)*N);
    pc.q_polys= (MLE_POLY*)malloc(sizeof(MLE_POLY)*N);
    for(size_t i=0;i<N;i++){
        pc.p_polys[i]= p_arr[i]; // we assume caller won't free them? or we do a deep copy?
        pc.q_polys[i]= q_arr[i];
    }
    pc.claim= claim;
    pc.num_vars= num_vars;
    pc.has_cached=0;
    return pc;
}

void prodcheck_free(ProdCheck pc){
    if(pc.p_polys){
        // optionally mlp_free each 
        for(size_t i=0;i<pc.N;i++){ // free them if we made a deep copy
            // mlp_free(&pc.p_polys[i].coeffs);
            // mlp_free(&pc.q_polys[i]);
        }
        free(pc.p_polys);
        free(pc.q_polys);
    }
    points_free(pc.challenges);
    pc.p_polys=NULL;
    pc.q_polys=NULL;
    pc.N=0;
    pc.has_cached=0;
}


CompressedPoly* prodcheck_round_polynomial(ProdCheck* pc)
{
    size_t p0_len= pc->p_polys[0].len;
    assert(p0_len>1 && "The protocol is already complete");

    // if has_cached => return cached_round_msg
    if(pc->has_cached){
        return pc->cached_round_msg;
    }

    size_t half= p0_len/2;

    // We'll do a local array of [3] for partial sums, then combine them
    F128 accum[3];
    accum[0]= f128_zero(); // sum of pq_zero across all i
    accum[1]= f128_zero(); // sum of pq_one
    accum[2]= f128_zero(); // sum of pq_inf

    for(size_t i=0; i< half;i++){
        F128 loc0= f128_zero(); 
        F128 loc1= f128_zero();
        F128 loc2= f128_zero();
        // for each j in [0..N)
        for(size_t j=0;j<pc->N;j++){
            F128 p_low= mlp_get(&pc->p_polys[j], 2*i);
            F128 p_high= mlp_get(&pc->p_polys[j], 2*i+1);
            F128 q_low= mlp_get(&pc->q_polys[j], 2*i);
            F128 q_high= mlp_get(&pc->q_polys[j], 2*i+1);

            // pq_zero += p_low * q_low
            loc0= f128_add(loc0, f128_mul(p_low, q_low));
            // pq_one += p_high * q_high
            loc1= f128_add(loc1, f128_mul(p_high, q_high));
            // pq_inf += (p_low+p_high)*(q_low+q_high)
            F128 sumP= f128_add(p_low, p_high);
            F128 sumQ= f128_add(q_low, q_high);
            loc2= f128_add(loc2, f128_mul(sumP, sumQ));
        }
        // reduce => accum
        accum[0]= f128_add(accum[0], loc0);
        accum[1]= f128_add(accum[1], loc1);
        accum[2]= f128_add(accum[2], loc2);
    }

    // poly[1] += poly[0] + poly[2];
    // define poly array:
    Points* tri;
    tri= (Points*)malloc(sizeof(Points));
    tri->len= 3;
    tri->elems= (F128*)malloc(sizeof(F128)*3);

    tri->elems[0]= accum[0];
    tri->elems[1]= accum[1];
    tri->elems[2]= accum[2];

    // poly[1] += poly[0] + poly[2]
    tri->elems[1]= f128_add(tri->elems[1], f128_add(tri->elems[0], tri->elems[2]));

    // compress => c0= tri[0], c1= tri[1], computed_claim= tri[2]
    CompressedPoly* cr = compress_poly(tri);
    free(tri->elems);
    free(tri);
    
    // ensure computed_claim == pc->claim
    if(!f128_eq(cr->sum, pc->claim)){
        fprintf(stderr,"Claim does not match expected value.\n");
        exit(1);
    }

    // store in cache
    pc->cached_round_msg= cr;
    pc->has_cached=1;

    return cr;
}

void prodcheck_bind(ProdCheck* pc, F128 r, int challenge_index)
{
    size_t p0_len= pc->p_polys[0].len;
    assert(p0_len>1 && "The protocol is already complete");

    size_t half= p0_len/2;

    // first retrieve the "round_poly" => do the same logic as above
    CompressedPoly* cp = prodcheck_round_polynomial(pc);
    UnivariatePolynomial* poly = uncompress_poly(cp);
    F128 new_claim = polynomial_evaluate_at(poly, r);

    pc->claim= new_claim;

    // push r into pc->challenges
    pc->challenges->elems[challenge_index]= r;

    // We now build new polynomials p_new, q_new each of length= half
    MLE_POLY* p_new = (MLE_POLY*)malloc(sizeof(MLE_POLY)* pc->N);
    MLE_POLY* q_new = (MLE_POLY*)malloc(sizeof(MLE_POLY)* pc->N);

    for(size_t i=0;i< pc->N;i++){
        p_new[i].coeffs= (F128*)malloc(sizeof(F128)* half);
        p_new[i].len= half;
        q_new[i].coeffs= (F128*)malloc(sizeof(F128)* half);
        q_new[i].len= half;
    }

    // halving
    for(size_t i=0;i< pc->N;i++){
        const MLE_POLY *p= &pc->p_polys[i];
        const MLE_POLY *q= &pc->q_polys[i];
        for(size_t j=0;j< half;j++){
            F128 p_low= mlp_get(p, 2*j);
            F128 p_high= mlp_get(p, 2*j+1);
            F128 q_low= mlp_get(q, 2*j);
            F128 q_high= mlp_get(q, 2*j+1);

            // p_new[i][j]= p_low + (p_high+p_low)* r
            F128 sumP= f128_add(p_high, p_low);
            F128 pr= f128_mul(sumP, r);
            F128 p_val= f128_add(p_low, pr);

            p_new[i].coeffs[j]= p_val;

            // q_new
            F128 sumQ= f128_add(q_high, q_low);
            F128 qr= f128_mul(sumQ, r);
            F128 q_val= f128_add(q_low, qr);
            q_new[i].coeffs[j]= q_val;
        }
    }

    // free old polynomials
    // TODO BETTER CLEAN UP HERE
    // for(size_t i=0; i< pc->N;i++){
    //     mlp_free(&pc->p_polys[i]); 
    //     mlp_free(&pc->q_polys[i]);
    // }
    free(pc->p_polys);
    free(pc->q_polys);

    pc->p_polys= p_new;
    pc->q_polys= q_new;

    pc->has_cached=0; // clear the cached msg
    pc->cached_round_msg= NULL; // TODO should we free this?
}

// "finish(self) -> Output"
ProdCheckOutput prodcheck_finish(ProdCheck pc)
{
    // build p_evaluations => array of length pc.N
    F128* pvals = (F128*)malloc(sizeof(F128)* pc.N);
    F128* qvals = (F128*)malloc(sizeof(F128)* pc.N);

    for(size_t i=0; i< pc.N; i++){
        // ensure length==1
        if(pc.p_polys[i].len !=1){
            fprintf(stderr,"The protocol is not complete\n");
            exit(1);
        }
        if(pc.q_polys[i].len !=1){
            fprintf(stderr,"The protocol is not complete\n");
            exit(1);
        }
        pvals[i]= pc.p_polys[i].coeffs[0];
        qvals[i]= pc.q_polys[i].coeffs[0];
    }

    Evaluations* fe_p= malloc(sizeof(Evaluations));
    Evaluations* fe_q= malloc(sizeof(Evaluations));
    fe_p->len= pc.N;
    fe_q->len= pc.N;
    fe_p->elems= pvals;
    fe_q->elems= qvals;

    // free(pvals);
    // free(qvals);

    ProdCheckOutput out;
    out.p_evaluations= fe_p;
    out.q_evaluations= fe_q;
    return out;
}