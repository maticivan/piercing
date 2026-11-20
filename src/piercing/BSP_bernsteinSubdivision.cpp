//*************************************************************************************************
//*************************************************************************************************
//* The MIT License (MIT)                                                                         *
//* Copyright (C) 2026                                                                            *
//* Ivan Matic and Rados Radoicic                                                                 *
//*                                                                                               *
//* Permission is hereby granted, free of charge, to any person obtaining a copy of this          *
//* software and associated documentation files (the "Software"), to deal in the Software         *
//* without restriction, including without limitation the rights to use, copy, modify, merge,     *
//* publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons    *
//* to whom the Software is furnished to do so, subject to the following conditions:              *
//*                                                                                               *
//* The above copyright notice and this permission notice shall be included in all copies or      *
//* substantial portions of the Software.                                                         *
//*                                                                                               *
//* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,           *
//* INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR      *
//* PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE     *
//* FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR          *
//* OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        *
//* DEALINGS IN THE SOFTWARE.                                                                     *
//*************************************************************************************************

// Bernstein expansions and the subdivision proof of non-negativity on [0,1]^2
// (Section 4 of the paper).
//
// A polynomial g of degree at most n in z1 and at most m in z2 is stored as an
// (n+1) x (m+1) array of integers. Two kinds of arrays occur:
//   - monomial coefficients:  g = sum_{k,l} p[k][l] z1^k z2^l,
//   - Bernstein coefficients: g = sum_{i,j} c[i][j] B^n_i(z1) B^m_j(z2).
// Multiplying all coefficients by the same positive integer does not change the sign of
// g at any point, nor the sign of any coefficient. The functions below use this freely,
// so that all arithmetic is in integers.
//
// All integers are fixed-width and "checked": an operation whose exact result does not
// fit throws std::overflow_error. The main program catches it and then refuses to
// print the verdict. No result is ever rounded.

#ifndef _INCL_BSP_CPP
#define _INCL_BSP_CPP

namespace mp = boost::multiprecision;
using Int1024 = mp::number<mp::cpp_int_backend<1024, 1024,
                    mp::signed_magnitude, mp::checked, void>>;

namespace BSP{

template<typename BigNum>
using Matrix=std::vector<std::vector<BigNum> >;

template<typename BigNum>
Matrix<BigNum> zeroMatrix(long rows, long cols){
    Matrix<BigNum> res;
    res.resize(rows);
    for(long i=0;i<rows;++i){
        res[i].resize(cols);
        for(long j=0;j<cols;++j){res[i][j]=0;}
    }
    return res;
}

template<typename BigNum>
BigNum absValue(const BigNum& x){
    if(x<0){return -x;}
    return x;
}

template<typename BigNum>
BigNum gcd(BigNum a, BigNum b){
    a=absValue(a);b=absValue(b);
    while(b!=0){
        BigNum r=a%b;
        a=b;b=r;
    }
    return a;
}

template<typename BigNum>
BigNum lcm(const BigNum& a, const BigNum& b){
    if((a==0)||(b==0)){return 0;}
    return absValue(a/gcd(a,b)*b);
}

// binomial[n][k] = C(n,k) for 0<=k<=n<=maxN, by Pascal's rule.
template<typename BigNum>
Matrix<BigNum> binomialTable(long maxN){
    Matrix<BigNum> res=zeroMatrix<BigNum>(maxN+1,maxN+1);
    for(long n=0;n<=maxN;++n){
        res[n][0]=1;
        for(long k=1;k<=n;++k){
            res[n][k]=res[n-1][k-1]+res[n-1][k];
        }
    }
    return res;
}

// Converts monomial coefficients to Bernstein coefficients of the same bidegree (n,m),
// multiplied by the positive integer L_n*L_m, where L_n = lcm of C(n,0),...,C(n,n).
// It uses   z^k = sum_{i=k}^{n} [C(i,k)/C(n,k)] B^n_i(z),
// so that   c[i][j] = sum_{k<=i, l<=j} C(i,k) C(j,l) / (C(n,k) C(m,l)) p[k][l].
template<typename BigNum>
Matrix<BigNum> monomialToBernsteinTimesPositiveInteger(const Matrix<BigNum>& p){
    long n=static_cast<long>(p.size())-1;
    long m=static_cast<long>(p[0].size())-1;
    long maxN=(n>m)?n:m;
    Matrix<BigNum> C=binomialTable<BigNum>(maxN);
    BigNum Ln=1,Lm=1;
    for(long k=0;k<=n;++k){Ln=lcm(Ln,C[n][k]);}
    for(long l=0;l<=m;++l){Lm=lcm(Lm,C[m][l]);}
    // q[k][l] = p[k][l] * (Ln/C(n,k)) * (Lm/C(m,l)); the divisions are exact.
    Matrix<BigNum> q=zeroMatrix<BigNum>(n+1,m+1);
    for(long k=0;k<=n;++k){
        for(long l=0;l<=m;++l){
            q[k][l]=p[k][l]*(Ln/C[n][k])*(Lm/C[m][l]);
        }
    }
    Matrix<BigNum> c=zeroMatrix<BigNum>(n+1,m+1);
    for(long i=0;i<=n;++i){
        for(long j=0;j<=m;++j){
            BigNum s=0;
            for(long k=0;k<=i;++k){
                for(long l=0;l<=j;++l){
                    if(q[k][l]!=0){s+=C[i][k]*C[j][l]*q[k][l];}
                }
            }
            c[i][j]=s;
        }
    }
    return c;
}

// One variable. Given Bernstein coefficients c[0..n] of a polynomial on [0,1], returns the
// Bernstein coefficients on [0,1/2] and on [1/2,1], each rescaled to [0,1] and multiplied
// by 2^n. This is de Casteljau's algorithm at 1/2, with sums in place of averages:
// after r steps the first and last entries are 2^r times the r-th left coefficient and
// the (n-r)-th right coefficient.
template<typename BigNum>
std::pair<std::vector<BigNum>,std::vector<BigNum> > splitAtOneHalf(const std::vector<BigNum>& c){
    long n=static_cast<long>(c.size())-1;
    std::vector<BigNum> t=c;
    std::vector<BigNum> left,right;
    left.resize(n+1);right.resize(n+1);
    BigNum twoPower=1;
    for(long r=0;r<n;++r){twoPower*=2;}   // 2^n
    left[0]=t[0]*twoPower;
    right[n]=t[n]*twoPower;
    for(long r=1;r<=n;++r){
        for(long i=0;i+r<=n;++i){t[i]=t[i]+t[i+1];}
        twoPower/=2;                        // 2^(n-r)
        left[r]=t[0]*twoPower;
        right[n-r]=t[n-r]*twoPower;
    }
    return {left,right};
}

// Two variables, split in z1. Returns the Bernstein coefficients of g on
// [0,1/2]x[0,1] and on [1/2,1]x[0,1], each rescaled to [0,1]^2 and multiplied by 2^n.
// Each column of c is the list of coefficients in z1 for one fixed index in z2.
template<typename BigNum>
std::pair<Matrix<BigNum>,Matrix<BigNum> > splitColumnsAtOneHalf(const Matrix<BigNum>& c){
    long n=static_cast<long>(c.size())-1;
    long m=static_cast<long>(c[0].size())-1;
    Matrix<BigNum> low=zeroMatrix<BigNum>(n+1,m+1), high=zeroMatrix<BigNum>(n+1,m+1);
    for(long j=0;j<=m;++j){
        std::vector<BigNum> column;column.resize(n+1);
        for(long i=0;i<=n;++i){column[i]=c[i][j];}
        std::pair<std::vector<BigNum>,std::vector<BigNum> > lr=splitAtOneHalf(column);
        for(long i=0;i<=n;++i){low[i][j]=lr.first[i];high[i][j]=lr.second[i];}
    }
    return {low,high};
}

// Two variables, split in z2. Returns the Bernstein coefficients of g on
// [0,1]x[0,1/2] and on [0,1]x[1/2,1], each rescaled to [0,1]^2 and multiplied by 2^m.
// Each row of c is the list of coefficients in z2 for one fixed index in z1.
template<typename BigNum>
std::pair<Matrix<BigNum>,Matrix<BigNum> > splitRowsAtOneHalf(const Matrix<BigNum>& c){
    long n=static_cast<long>(c.size())-1;
    Matrix<BigNum> low,high;
    low.resize(n+1);high.resize(n+1);
    for(long i=0;i<=n;++i){
        std::pair<std::vector<BigNum>,std::vector<BigNum> > lr=splitAtOneHalf(c[i]);
        low[i]=lr.first;high[i]=lr.second;
    }
    return {low,high};
}

// Two variables. Returns the Bernstein coefficients of g on the four squares of half the
// side, each rescaled to [0,1]^2 and multiplied by a positive integer. The order is
// [0,1/2]x[0,1/2], [0,1/2]x[1/2,1], [1/2,1]x[0,1/2], [1/2,1]x[1/2,1].
template<typename BigNum>
std::vector<Matrix<BigNum> > splitIntoFour(const Matrix<BigNum>& c){
    std::pair<Matrix<BigNum>,Matrix<BigNum> > z1Halves=splitColumnsAtOneHalf(c);
    std::pair<Matrix<BigNum>,Matrix<BigNum> > lowQuarters=splitRowsAtOneHalf(z1Halves.first);
    std::pair<Matrix<BigNum>,Matrix<BigNum> > highQuarters=splitRowsAtOneHalf(z1Halves.second);
    return {lowQuarters.first,lowQuarters.second,highQuarters.first,highQuarters.second};
}

// Divides all coefficients by their (positive) greatest common divisor. This keeps the
// integers small and changes no sign.
template<typename BigNum>
void divideByContent(Matrix<BigNum>& c){
    BigNum g=0;
    for(long i=0;i<static_cast<long>(c.size());++i){
        for(long j=0;j<static_cast<long>(c[i].size());++j){g=gcd(g,c[i][j]);}
    }
    if(g>1){
        for(long i=0;i<static_cast<long>(c.size());++i){
            for(long j=0;j<static_cast<long>(c[i].size());++j){c[i][j]/=g;}
        }
    }
}

template<typename BigNum>
int allNonNegative(const Matrix<BigNum>& c){
    for(long i=0;i<static_cast<long>(c.size());++i){
        for(long j=0;j<static_cast<long>(c[i].size());++j){
            if(c[i][j]<0){return 0;}
        }
    }
    return 1;
}

struct SubdivisionResult{
public:
    int proved;              // 1 if g >= 0 on [0,1]^2 was proved
    long finalSquares;       // squares in the final list (all coefficients >= 0)
    long squaresExamined;    // all squares whose coefficients were computed
    long maxDepth;           // largest number of halvings of a side
    std::string failure;     // reason, if proved==0
    SubdivisionResult(){proved=0;finalSquares=0;squaresExamined=0;maxDepth=0;}
};

// Proves g >= 0 on [0,1]^2 from its Bernstein coefficients c (Section 4).
// A square on which all coefficients are >= 0 goes to the final list. A square with a
// negative coefficient is divided into four. The procedure stops with failure if a corner
// coefficient is negative (it equals the value of g at that corner, so g < 0 there),
// or if a square has been halved more than depthLimit times.
template<typename BigNum>
SubdivisionResult proveNonNegativeOnUnitSquare(const Matrix<BigNum>& c, long depthLimit){
    SubdivisionResult res;
    fst::stack<std::pair<Matrix<BigNum>,long> > toExamine;
    toExamine.push({c,0});
    while(!toExamine.empty()){
        std::pair<Matrix<BigNum>,long> current=toExamine.top();
        toExamine.pop();
        const Matrix<BigNum>& d=current.first;
        long depth=current.second;
        ++res.squaresExamined;
        if(depth>res.maxDepth){res.maxDepth=depth;}
        if(allNonNegative(d)){
            ++res.finalSquares;
            continue;
        }
        long n=static_cast<long>(d.size())-1, m=static_cast<long>(d[0].size())-1;
        if((d[0][0]<0)||(d[0][m]<0)||(d[n][0]<0)||(d[n][m]<0)){
            res.failure="the polynomial is negative at a corner of a square";
            return res;
        }
        if(depth>=depthLimit){
            res.failure="depth limit reached";
            return res;
        }
        std::vector<Matrix<BigNum> > children=splitIntoFour(d);
        for(long q=0;q<4;++q){
            divideByContent(children[q]);
            toExamine.push({children[q],depth+1});
        }
    }
    res.proved=1;
    return res;
}

} 
#endif
