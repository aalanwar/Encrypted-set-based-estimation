/*
 * Copyright 2013-2015 Raphael Bost
 *
 * This file is part of ciphermed.

 *  ciphermed is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ciphermed is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ciphermed.  If not, see <http://www.gnu.org/licenses/>. 2
 *
 */

#include "num_th_alg.hh"
#include "math_util.hh"
#include "util_gmp_rand.h"
#include "prime_seq.hh"
#include "mpz_class.hh"
#include <NTL/ZZ.h>
#include <cassert>

#include <iostream>

/* The algorithms here are from Victor Shoup's Book
 * A Computational Introduction to Number Theory and Algebra
 * or Shoup's NTL library (for Sophie Germain primes)
 */

std::vector<mpz_class> gen_rand_non_increasing_seq(const mpz_class &m, gmp_randstate_t state)
{
    std::vector<mpz_class> seq;

    mpz_class n = m;
    mpz_class new_n;
    do {
        // pick a new n between 1...n
        mpz_urandomm(new_n.get_mpz_t(),state,n.get_mpz_t());
        n = new_n + 1;
        seq.push_back(n);
    } while (n != 1);

    return seq;
}

std::vector<mpz_class> extract_prime_seq(const std::vector<mpz_class> &seq, int reps)
{
    std::vector<mpz_class> primes;
    for (size_t i = 0; i < seq.size(); i++) {
        if (mpz_class_probab_prime_p(seq[i],reps) != 0) {
            primes.push_back(seq[i]);
        }
    }

    return primes;
}

// generates a random integer with its factorization (returns the factorization)
std::vector<mpz_class> gen_rand_number_factorization(const mpz_class &m, mpz_class *result, gmp_randstate_t state, int reps)
{
    for (; ; ) {
        std::vector<mpz_class> seq = gen_rand_non_increasing_seq(m,state);
        std::vector<mpz_class> primes = extract_prime_seq(seq, reps);

        mpz_class y = 1;
        for (size_t i = 0; i < primes.size(); i++) {
            y *= primes[i];
            if (y > m) {
                break;
            }
        }
        if (y > m) {
            continue;
        }

        mpz_class x;
        mpz_urandomm(x.get_mpz_t(),state,m.get_mpz_t());
        x += 1;

        if (x <= y) {
            if (result) {
                *result = y;
            }
            return primes;
        }
    }
}

// generates a random prime p with the factorization of p-1 (returns the factorization)
std::vector<mpz_class> gen_rand_prime_with_factorization(const mpz_class &m, mpz_class *p, gmp_randstate_t state, int reps)
{
    for (; ; ) {
        mpz_class n;
        std::vector<mpz_class> factorization = gen_rand_number_factorization(m,&n,state,reps);

        if (mpz_class_probab_prime_p(n+1, reps) != 0) {
            if (p) {
                *p = n+1;
            }
            return factorization;
        }
    }
}

mpz_class simple_safe_prime_gen(size_t n_bits, gmp_randstate_t state, int reps)
{
    for (size_t count = 1; ; count ++) {
        mpz_class n;
        mpz_urandom_len(n.get_mpz_t(),state,n_bits);

        if (mpz_class_probab_prime_p(n,reps) !=0) {
            if (mpz_class_probab_prime_p(2*n+1,reps) != 0 ) {
                std::cout << count << " iterations needed to generate safe prime" << std::endl;
                return n;
            }
        }
    }
}


static long bit_count(long a)
{
    unsigned long aa;
    if (a < 0)
    aa = - ((unsigned long) a);
    else
    aa = a;

    long k = 0;
    while (aa) {
        k++;
        aa = aa >> 1;
    }

    return k;
}

/* The following code is just NTL's code for generating Germain primes using GMP */

// prime_bound computes a reasonable bound for trial
// division in the Miller-Rabin test.
// It is computed a bit on the "low" side, since being a bit
// low doesn't hurt much, but being too high can hurt a lot.

static
long prime_bound(long bn)
{
    long wn = (bn+NBITS_MAX-1)/NBITS_MAX;

    long fn;

    if (wn <= 36)
    fn = wn/4 + 1;
    else
    fn = long(1.67*sqrt(double(wn)));

    long prime_bnd;

    if (bit_count(bn) + bit_count(fn) > NBITS_MAX)
        prime_bnd = (1L << NBITS_MAX);
    else
        prime_bnd = bn*fn;

    return prime_bnd;
}

static
long ErrBoundTest(long kk, long tt, long nn)

{
    const double fudge = (1.0 + 1024.0/NTL_FDOUBLE_PRECISION);
    const double log2_3 = log2(3.0);
    const double log2_7 = log2(7.0);
    const double log2_20 = log2(20.0);

    double k = kk;
    double t = tt;
    double n = nn;

    if (k < 3 || t < 1) return 0;
    if (n < 1) return 1;

    // the following test is largely academic
    assert(9*t < NTL_FDOUBLE_PRECISION);

    double log2_k = log2(k);

    if ((n + log2_k)*fudge <= 2*t)
    return 1;

    if ((2*log2_k + 4.0 + n)*fudge <= 2*sqrt(k))
    return 2;

    if ((t == 2 && k >= 88) || (3 <= t && 9*t <= k && k >= 21)) {
        if ((1.5*log2_k + t + 4.0 + n)*fudge <= 0.5*log2(t) + 2*(sqrt(t*k)))
        return 3;
    }

    if (k <= 9*t && 4*t <= k && k >= 21) {
        if ( ((log2_3 + log2_7 + log2_k + n)*fudge <= log2_20 + 5*t)  &&
            ((log2_3 + (15.0/4.0)*log2_k + n)*fudge <= log2_7 + k/2 + 2*t) &&
            ((2*log2_3 + 2 + log2_k + n)*fudge <= k/4 + 3*t) )
        return 4;
    }

    if (4*t >= k && k >= 21) {
        if (((15.0/4.0)*log2_k + n)*fudge <= log2_7 + k/2 + 2*t)
        return 5;
    }

    return 0;
}

static long make_odd(mpz_class &n)
{
    long k = 0;

    while (mpz_even_p(n.get_mpz_t())) {
        n >>= 1;
        k++;
    }

    return k;
}

static long is_Miller_witness(const mpz_class& n, const mpz_class& x)
{
    mpz_class m(0), y(0), z(0);
    long j, k;

    if (x == 0) return 0;

    m = n-1;
    k = make_odd(m);

    z = mpz_class_powm(x,m,n);

    if (z == 1) return 0;

    j = 0;
    do {
        y = z;
        z = (y*y) %n;
        j++;
    } while (j != k && z != 1);

    if (z != 1) return 1;
    y = y + 1;
    if (y != n) return 1;
    return 0;
}

void gen_germain_prime(mpz_class& n, long k, gmp_randstate_t state, long err)
{
    assert(k > 1);
    assert(k <= (1L << 20));

    if (err < 1) err = 1;
    if (err > 512) err = 512;

    if (k == 2) {
        if (gmp_urandomm_ui(state,2))
        n = 3;
        else
        n = 2;

        return;
    }


    long prime_bnd = prime_bound(k);

    if (bit_count(prime_bnd) >= k/2)
    prime_bnd = (1L << (k/2-1));


    mpz_class two;
    two = 2;

    mpz_class n1;


    PrimeSeq s;

    mpz_class iter;
    iter = 0;


    for (;;) {
        iter++;

        mpz_urandom_len(n.get_mpz_t(),state,k);

        if (mpz_even_p(n.get_mpz_t())) {
            n = n+1;
        }

        s.reset(3);
        long p;

        long sieve_passed = 1;

        p = s.next();
        while (p && p < prime_bnd) {
            mpz_class r;
            mpz_tdiv_r_ui(r.get_mpz_t(),n.get_mpz_t(),p);

            if (r == 0) {
                sieve_passed = 0;
                break;
            }

            // test if 2*r + 1 = 0 (mod p)
            if (r == p-r-1) {
                sieve_passed = 0;
                break;
            }

            p = s.next();
        }

        if (!sieve_passed) continue;


        if (is_Miller_witness(n, two)) continue;

        n1 = 2*n+1;

        if (is_Miller_witness(n1, two)) continue;

        // now do t M-R iterations...just to make sure

        // First compute the appropriate number of M-R iterations, t
        // The following computes t such that
        //       p(k,t)*8/k <= 2^{-err}/(5*iter^{1.25})
        // which suffices to get an overall error probability of 2^{-err}.
        // Note that this method has the advantage of not requiring
        // any assumptions on the density of Germain primes.
        long iter_n_bits = mpz_sizeinbase(iter.get_mpz_t(),2);
        long err1 = std::max(1L, err + 7 + (5*iter_n_bits + 3)/4 - bit_count(k));
        long t;
        t = 1;
        while (!ErrBoundTest(k, t, err1))
        t++;

        if(mpz_probab_prime_p(n.get_mpz_t(),t))
            break;
    }
}

// Constructs a generator for the cyclic group \Z^*_p where p is a Sophie Germain prime
mpz_class get_generator_for_cyclic_group(const mpz_class &p, gmp_randstate_t state)
{
    mpz_class q = (p >> 1);
    mpz_class g;


    // find a generator for ZZ*_p
    // Shoup's algorithm

    mpz_class alpha, beta;
    do {
        mpz_urandomm(alpha.get_mpz_t(),state,p.get_mpz_t());

        beta = mpz_class_powm(alpha,q,p);
    } while (beta == 1);

    g = beta;

    do {
        mpz_urandomm(alpha.get_mpz_t(),state,p.get_mpz_t());
        beta = (alpha*alpha) %p;

    } while (beta == 1);

    g = (g*beta) %p;

    return g;
}
