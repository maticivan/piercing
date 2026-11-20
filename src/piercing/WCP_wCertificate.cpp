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

// The polynomial W_14 and the three statements (V1), (V2), (V3) of Section 4.
//
// The input file gives integers bNum[i][j] for 0<=i<=j<=14 and a positive integer bDen.
// The program sets bNum[j][i]=bNum[i][j], so that b_ij = bNum[i][j]/bDen = b_ji, and
//     W_14(x,y) = sum_{i,j} b_ij B^14_i(S) B^14_j(U),   S=(2x-1)^2, U=(2y-1)^2.     (3.1)
// Throughout, the program works with bDen*W_14, which has integer coefficients.
//
// (V1)  W_14 >= 0 on [0,1]^2.
//       As a polynomial in (S,U), W_14 has Bernstein coefficients b_ij, and (x,y)->(S,U)
//       maps [0,1]^2 onto itself. So (V1) is the subdivision proof for the array bNum.
// (V2)  The integral of W_14 over [0,1]^2 equals massNum/massDen.
//       W_14 = sum w_kl S^k U^l and the integral of S^k U^l is 1/((2k+1)(2l+1)).
// (V3)  beta - F_A >= 0 and beta - F_B >= 0 on [0,1]^2, beta = betaNum/betaDen, where
//       F_A(a,b) = (1+b-a) int_0^1 W_14(tau, (1-tau)a + tau b) dtau,          (3.2)
//       F_B(a,b) = (1-a+b) int_0^1 W_14(b tau, (1-tau)a + tau) dtau.          (3.3)

#ifndef _INCL_WCP_CPP
#define _INCL_WCP_CPP

namespace WCP{

const long DEG=14;   // W_14 has degree 14 in S and in U

// ------------------------------------------------------------------------------ input

// Parses a decimal integer, with an optional sign. Anything else is an error.
template<typename BigNum>
std::pair<BigNum,int> parseInteger(const std::string& s){
    std::pair<BigNum,int> res;res.first=0;res.second=0;
    long pos=0,len=static_cast<long>(s.size());
    while((pos<len)&&std::isspace(static_cast<unsigned char>(s[pos]))){++pos;}
    while((len>pos)&&std::isspace(static_cast<unsigned char>(s[len-1]))){--len;}
    int negative=0;
    if((pos<len)&&((s[pos]=='-')||(s[pos]=='+'))){negative=(s[pos]=='-');++pos;}
    if(pos>=len){return res;}
    BigNum value=0;
    for(;pos<len;++pos){
        if((s[pos]<'0')||(s[pos]>'9')){return res;}
        value=value*10+(s[pos]-'0');
    }
    res.first=negative?(-value):value;
    res.second=1;
    return res;
}

// Reads a value block {{n1}{n2}...}{{...}}... into rows of integers.
template<typename BigNum>
std::pair<std::vector<std::vector<BigNum> >,int> readRows(const std::string& block){
    std::pair<std::vector<std::vector<BigNum> >,int> res;res.second=1;
    std::vector<std::string> rows=SF::stringToVector(block,"{","}");
    res.first.resize(rows.size());
    for(long r=0;r<static_cast<long>(rows.size());++r){
        std::vector<std::string> entries=SF::stringToVector(rows[r],"{","}");
        for(long e=0;e<static_cast<long>(entries.size());++e){
            std::pair<BigNum,int> x=parseInteger<BigNum>(entries[e]);
            if(x.second==0){res.second=0;}
            res.first[r].push_back(x.first);
        }
    }
    return res;
}

template<typename BigNum>
struct InputData{
public:
    BSP::Matrix<BigNum> bNum;   // 15 x 15, symmetric
    BigNum bDen, betaNum, betaDen, massNum, massDen;
    int ok;
    std::string error;
};

template<typename BigNum>
InputData<BigNum> readInput(const std::string& fileContent){
    InputData<BigNum> d;d.ok=0;
    std::map<std::string,std::string> records=SF::stringToMap(fileContent,"[key]","[/key]","[value]","[/value]");
    std::vector<std::string> names={"bUpperNum","bDen","betaNum","betaDen","massNum","massDen"};
    for(long i=0;i<static_cast<long>(names.size());++i){
        if(records.find(names[i])==records.end()){d.error="missing record "+names[i];return d;}
    }
    std::vector<BigNum*> scalars={&d.bDen,&d.betaNum,&d.betaDen,&d.massNum,&d.massDen};
    for(long i=0;i<5;++i){
        std::pair<std::vector<std::vector<BigNum> >,int> v=readRows<BigNum>(records[names[i+1]]);
        if((v.second==0)||(v.first.size()!=1)||(v.first[0].size()!=1)){
            d.error="record "+names[i+1]+" must hold one integer";return d;
        }
        *(scalars[i])=v.first[0][0];
    }
    if((d.bDen<=0)||(d.betaNum<=0)||(d.betaDen<=0)||(d.massDen<=0)){
        d.error="bDen, betaNum, betaDen and massDen must be positive";return d;
    }
    std::pair<std::vector<std::vector<BigNum> >,int> rows=readRows<BigNum>(records["bUpperNum"]);
    if(rows.second==0){d.error="bUpperNum contains an entry that is not an integer";return d;}
    if(static_cast<long>(rows.first.size())!=DEG+1){d.error="bUpperNum must have 15 rows";return d;}
    d.bNum=BSP::zeroMatrix<BigNum>(DEG+1,DEG+1);
    for(long i=0;i<=DEG;++i){
        if(static_cast<long>(rows.first[i].size())!=DEG+1-i){
            d.error="row "+std::to_string(i)+" of bUpperNum must have "+std::to_string(DEG+1-i)+" entries";
            return d;
        }
        for(long j=i;j<=DEG;++j){
            d.bNum[i][j]=rows.first[i][j-i];
            d.bNum[j][i]=rows.first[i][j-i];      // b_ji = b_ij
        }
    }
    d.ok=1;
    return d;
}

// ------------------------------------------------------------------- W_14 in monomials

// Returns w with bDen*W_14 = sum_{k,l} w[k][l] S^k U^l. It uses
//   B^14_i(Z) = sum_{k=i}^{14} C(14,i) C(14-i,k-i) (-1)^(k-i) Z^k.
template<typename BigNum>
BSP::Matrix<BigNum> monomialFormOfW(const BSP::Matrix<BigNum>& bNum){
    BSP::Matrix<BigNum> C=BSP::binomialTable<BigNum>(DEG);
    BSP::Matrix<BigNum> beta=BSP::zeroMatrix<BigNum>(DEG+1,DEG+1);   // B_i = sum_k beta[i][k] Z^k
    for(long i=0;i<=DEG;++i){
        for(long k=i;k<=DEG;++k){
            beta[i][k]=C[DEG][i]*C[DEG-i][k-i];
            if((k-i)%2==1){beta[i][k]=-beta[i][k];}
        }
    }
    BSP::Matrix<BigNum> w=BSP::zeroMatrix<BigNum>(DEG+1,DEG+1);
    for(long i=0;i<=DEG;++i){
        for(long j=0;j<=DEG;++j){
            if(bNum[i][j]==0){continue;}
            for(long k=i;k<=DEG;++k){
                for(long l=j;l<=DEG;++l){
                    w[k][l]+=bNum[i][j]*beta[i][k]*beta[j][l];
                }
            }
        }
    }
    return w;
}

// ----------------------------------------------- polynomials in three variables (a,b,tau)

// Dense array of coefficients of sum c(i,j,k) a^i b^j tau^k, 0<=i<=da, 0<=j<=db, 0<=k<=dt.
template<typename BigNum>
struct Poly3{
public:
    long da,db,dt;
    std::vector<BigNum> c;
    Poly3(long _da=0,long _db=0,long _dt=0){
        da=_da;db=_db;dt=_dt;
        c.resize((da+1)*(db+1)*(dt+1));
        for(long i=0;i<static_cast<long>(c.size());++i){c[i]=0;}
    }
    BigNum& at(long i,long j,long k){return c[(i*(db+1)+j)*(dt+1)+k];}
    const BigNum& at(long i,long j,long k) const{return c[(i*(db+1)+j)*(dt+1)+k];}
};

template<typename BigNum>
Poly3<BigNum> constantPoly3(const BigNum& x){
    Poly3<BigNum> p(0,0,0);p.at(0,0,0)=x;return p;
}

template<typename BigNum>
Poly3<BigNum> operator*(const Poly3<BigNum>& p, const Poly3<BigNum>& q){
    Poly3<BigNum> r(p.da+q.da,p.db+q.db,p.dt+q.dt);
    for(long i=0;i<=q.da;++i){for(long j=0;j<=q.db;++j){for(long k=0;k<=q.dt;++k){
        const BigNum& y=q.at(i,j,k);
        if(y==0){continue;}
        for(long i1=0;i1<=p.da;++i1){for(long j1=0;j1<=p.db;++j1){for(long k1=0;k1<=p.dt;++k1){
            const BigNum& x=p.at(i1,j1,k1);
            if(x!=0){r.at(i+i1,j+j1,k+k1)+=x*y;}
        }}}
    }}}
    return r;
}

template<typename BigNum>
Poly3<BigNum> operator+(const Poly3<BigNum>& p, const Poly3<BigNum>& q){
    Poly3<BigNum> r(std::max(p.da,q.da),std::max(p.db,q.db),std::max(p.dt,q.dt));
    for(long i=0;i<=p.da;++i){for(long j=0;j<=p.db;++j){for(long k=0;k<=p.dt;++k){r.at(i,j,k)+=p.at(i,j,k);}}}
    for(long i=0;i<=q.da;++i){for(long j=0;j<=q.db;++j){for(long k=0;k<=q.dt;++k){r.at(i,j,k)+=q.at(i,j,k);}}}
    return r;
}

// Returns bDen*W_14(X,Y) as a polynomial in (a,b,tau), for polynomials X, Y in (a,b,tau).
// With S=(2X-1)^2 and U=(2Y-1)^2:
//     bDen*W_14 = sum_k S^k Q_k,   Q_k = sum_l w[k][l] U^l,
// and the outer sum is evaluated by Horner's rule in S.
template<typename BigNum>
Poly3<BigNum> composeW(const BSP::Matrix<BigNum>& w, const Poly3<BigNum>& X, const Poly3<BigNum>& Y){
    Poly3<BigNum> minusOne=constantPoly3<BigNum>(-1), two=constantPoly3<BigNum>(2);
    Poly3<BigNum> twoXminus1=two*X+minusOne, twoYminus1=two*Y+minusOne;
    Poly3<BigNum> S=twoXminus1*twoXminus1, U=twoYminus1*twoYminus1;
    std::vector<Poly3<BigNum> > Upow;
    Upow.push_back(constantPoly3<BigNum>(1));
    for(long l=1;l<=DEG;++l){Upow.push_back(Upow[l-1]*U);}
    std::vector<Poly3<BigNum> > Q;
    for(long k=0;k<=DEG;++k){
        Poly3<BigNum> Qk;
        for(long l=0;l<=DEG;++l){
            if(w[k][l]!=0){Qk=Qk+constantPoly3<BigNum>(w[k][l])*Upow[l];}
        }
        Q.push_back(Qk);
    }
    Poly3<BigNum> result=Q[DEG];
    for(long k=DEG-1;k>=0;--k){result=result*S+Q[k];}
    return result;
}

// Given G(a,b,tau) with integer coefficients, returns the integer polynomial
//     L * (prefactor)(a,b) * int_0^1 G(a,b,tau) dtau,
// where L = lcm(1,...,dt+1), using int_0^1 tau^k dtau = 1/(k+1).
template<typename BigNum>
BSP::Matrix<BigNum> integrateTauTimesL(const Poly3<BigNum>& G, const Poly3<BigNum>& prefactor, BigNum& L){
    L=1;
    for(long k=0;k<=G.dt;++k){L=BSP::lcm(L,BigNum(k+1));}
    Poly3<BigNum> H(G.da,G.db,0);
    for(long i=0;i<=G.da;++i){for(long j=0;j<=G.db;++j){
        BigNum s=0;
        for(long k=0;k<=G.dt;++k){s+=G.at(i,j,k)*(L/(k+1));}
        H.at(i,j,0)=s;
    }}
    Poly3<BigNum> P=H*prefactor;
    BSP::Matrix<BigNum> res=BSP::zeroMatrix<BigNum>(P.da+1,P.db+1);
    for(long i=0;i<=P.da;++i){for(long j=0;j<=P.db;++j){res[i][j]=P.at(i,j,0);}}
    return res;
}

// Returns the monomial coefficients of the integer polynomial
//     betaDen*bDen*L*(beta - F)(a,b) = betaNum*bDen*L - betaDen*(bDen*L*F)(a,b),
// a positive multiple of beta-F, where F is (3.2) or (3.3).
template<typename BigNum>
BSP::Matrix<BigNum> betaMinusF(const InputData<BigNum>& d, const BSP::Matrix<BigNum>& w, char family){
    Poly3<BigNum> a(1,0,0),b(0,1,0),tau(0,0,1),one=constantPoly3<BigNum>(1);
    a.at(1,0,0)=1;b.at(0,1,0)=1;tau.at(0,0,1)=1;
    Poly3<BigNum> minusA=constantPoly3<BigNum>(-1)*a, minusTau=constantPoly3<BigNum>(-1)*tau;
    Poly3<BigNum> X,Y,prefactor;
    if(family=='A'){
        X=tau;                                    // (tau, (1-tau)a + tau b)
        Y=(one+minusTau)*a+tau*b;
        prefactor=one+b+minusA;                   // 1+b-a
    }
    else{
        X=b*tau;                                  // (b tau, (1-tau)a + tau)
        Y=(one+minusTau)*a+tau;
        prefactor=one+minusA+b;                   // 1-a+b
    }
    Poly3<BigNum> G=composeW(w,X,Y);              // bDen*W_14(X,Y)
    BigNum L;
    BSP::Matrix<BigNum> P=integrateTauTimesL(G,prefactor,L);   // bDen*L*F
    BSP::Matrix<BigNum> res=P;
    for(long i=0;i<static_cast<long>(res.size());++i){
        for(long j=0;j<static_cast<long>(res[i].size());++j){
            res[i][j]=-d.betaDen*P[i][j];
        }
    }
    res[0][0]+=d.betaNum*d.bDen*L;
    return res;
}

// ------------------------------------------------------------------- the three statements

template<typename BigNum>
void printSubdivision(const std::string& name, const BSP::SubdivisionResult& r){
    std::cout<<name<<": ";
    if(r.proved){
        std::cout<<"proved. Final list: "<<r.finalSquares<<" squares ("
                 <<r.squaresExamined<<" examined, at most "<<r.maxDepth<<" halvings).\n";
    }
    else{
        std::cout<<"NOT proved: "<<r.failure<<".\n";
    }
}

// (V2): sum w[k][l]/((2k+1)(2l+1)) / bDen == massNum/massDen, checked as
//     massDen * sum w[k][l]*(M/(2k+1))*(M/(2l+1)) == massNum * bDen * M^2,  M=lcm(1,3,...,29).
template<typename BigNum>
int checkMass(const InputData<BigNum>& d, const BSP::Matrix<BigNum>& w){
    BigNum M=1;
    for(long k=0;k<=DEG;++k){M=BSP::lcm(M,BigNum(2*k+1));}
    BigNum s=0;
    for(long k=0;k<=DEG;++k){for(long l=0;l<=DEG;++l){
        s+=w[k][l]*(M/(2*k+1))*(M/(2*l+1));
    }}
    return (d.massDen*s==d.massNum*d.bDen*M*M);
}

template<typename BigNum>
int mainWork(const std::string& fileContent){
    const long depthLimit=20;
    InputData<BigNum> d=readInput<BigNum>(fileContent);
    if(!d.ok){std::cout<<"Input error: "<<d.error<<"\n";return 0;}
    std::cout<<"Read the coefficients b_ij, 0<=i<=j<=14, and set b_ji=b_ij.\n";

    BSP::SubdivisionResult v1=BSP::proveNonNegativeOnUnitSquare(d.bNum,depthLimit);
    printSubdivision<BigNum>("(V1) W_14 >= 0 on [0,1]^2",v1);

    BSP::Matrix<BigNum> w=monomialFormOfW(d.bNum);
    int v2=checkMass(d,w);
    std::cout<<"(V2) integral of W_14 equals "<<d.massNum<<"/"<<d.massDen<<": "
             <<(v2?"proved.":"NOT proved.")<<"\n";

    BSP::Matrix<BigNum> pA=betaMinusF(d,w,'A');
    BSP::Matrix<BigNum> cA=BSP::monomialToBernsteinTimesPositiveInteger(pA);
    BSP::SubdivisionResult v3A=BSP::proveNonNegativeOnUnitSquare(cA,depthLimit);
    printSubdivision<BigNum>("(V3) beta - F_A >= 0 on [0,1]^2, bidegree ("
                             +std::to_string(pA.size()-1)+","+std::to_string(pA[0].size()-1)+")",v3A);

    BSP::Matrix<BigNum> pB=betaMinusF(d,w,'B');
    BSP::Matrix<BigNum> cB=BSP::monomialToBernsteinTimesPositiveInteger(pB);
    BSP::SubdivisionResult v3B=BSP::proveNonNegativeOnUnitSquare(cB,depthLimit);
    printSubdivision<BigNum>("(V3) beta - F_B >= 0 on [0,1]^2, bidegree ("
                             +std::to_string(pB.size()-1)+","+std::to_string(pB[0].size()-1)+")",v3B);

    int success=v1.proved&&v2&&v3A.proved&&v3B.proved;
    if(success){
        BigNum num=d.massNum*d.betaDen, den=d.massDen*d.betaNum;
        BigNum g=BSP::gcd(num,den);
        std::cout<<"alpha/beta = "<<(num/g)<<"/"<<(den/g)<<"\n";
    }
    return success;
}

} 
#endif
