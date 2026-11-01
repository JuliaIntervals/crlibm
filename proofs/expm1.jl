# expm1.c quick phase: common path (|x| >= 1/4, through exp) and direct path (|x| < 1/4)

module Expm1Proof
include("common.jl")

const H = "expm1.h"
const T1 = reshape(cdoubles(H, "twoPowerIndex1"), 3, 64); const T2 = reshape(cdoubles(H, "twoPowerIndex2"), 3, 64)
const CINV = cdouble(H, "log2InvMult2L")
const Lh, Lm = cdouble(H, "msLog2Div2Lh"), cdouble(H, "msLog2Div2Lm")
const C3, C4 = cdouble(H, "quickCommonpolyC3h"), cdouble(H, "quickCommonpolyC4h")
const D = [cdouble(H, "quickDirectpolyC$(k)h") for k in 3:9]
const OVF, MINUS1 = cdouble(H, "OVERFLOWBOUND"), cdouble(H, "MINUSONEBOUND")
@assert cint(H, "INDEXMASK1") == 0x3f && cint(H, "INDEXMASK2") == 0xfc0
const LN2 = log(Arb(2))

# |p - f| on discs |z| <= R: majorants of the difference series
expmaj(R, cs) = sum(abs(big(cs[k]) - 1/factorial(big(k)))*R^k for k in eachindex(cs)) +
                ub(exp(Arb(R))) - sum(R^k/factorial(big(k)) for k in 0:maximum(eachindex(cs)))

# ---------------------------------------------------------------- common path

const ETAB = let
    e1 = maximum(ub((Arb(T1[1, i+1]) + Arb(T1[2, i+1]))/Arb(2)^(Arb(i)/4096) - 1) for i in 0:63)
    e2 = maximum(ub((Arb(T2[1, i+1]) + Arb(T2[2, i+1]))/Arb(2)^(Arb(i)/64) - 1) for i in 0:63)
    (1 + e1)*(1 + e2)*(1 + EM22) - 1
end
const XMAX = max(abs(big(OVF)), abs(big(MINUS1)))
const KMAX = XMAX*CINV + 1
# rh + rm = x + Mul12(Lh, kd) + RN(RN(kd*Lm) + s2), x + s1 exact (Sterbenz)
const ERED = hu(KMAX*abs(Lm)) + hu(hu(KMAX*abs(Lh)) + KMAX*abs(Lm)*(1 + big(2.0)^-50)) + KMAX*ub(Arb(Lh) + Arb(Lm) + LN2/4096)
const RMAX = ((0.5 + hu(KMAX) + XMAX*ub(Arb(CINV) - 4096/LN2))*ub(LN2/4096) + ERED)*(1 + big(2.0)^-40)
@assert expmaj(1 + RMAX, Dict(3 => C3, 4 => C4)) <= 4

# error of exph + expm against exp(rh + rm) * (table value), in units of the table value
function commonR(R::Arb)
    Rm = ub(R); Hm = hu(Rm)
    rh = FQ(R)
    sq = fmul(rh, rh); c3 = fmul(FQ(C3), rh); sqh = fmul(FQ(0.5), sq; exact = true)
    hp = fadd(fmul(c3, sq), fmul(FQ(C4), fmul(sq, sq))); hpws = fadd(sqh, hp)
    t8e = hpws.e + hu(ub(hpws.v) + hpws.e + Hm)
    t9m = ub(R + hpws.v) + Hm + t8e
    t9e = t8e + hu(t9m)
    A = tm(s -> exp(s) - 1 - (s + s*s/2 + Arb(C3)*s*s*s + Arb(C4)*s*s*s*s), R, 4.0, 1.0; K = 14)
    rmt = Hm*Rm*1.001 + Hm^2
    ETAB*(1 + t9m) + ub(A) + rmt + t9e + 2*big(2.0)^-53*t9m + big(2.0)^-104 + ERED*1.0002
end
const EXPEPS = maximum(commonR(ival(a, b))/lb(exp(ival(a, b))) for (a, b) in linpieces(-RMAX, RMAX, 256))

function commonX(a::BigFloat, b::BigFloat)
    X = ival(a, b); E = exp(X)
    e = ub(E)*EXPEPS*(1 + big(2.0)^-50) + hu(hu(ub(E) + 1)*2 + hu(ub(E)))*2
    metrics(lb(expm1(X)), e)
end

# ---------------------------------------------------------------- direct path

# relative error of polyh + polyl for x in X; long: xIntHi > SPECIALINTERVALBOUND
function directpoly(X::Arb, long::Bool)
    xq = FQ(X); W = X*X; Wm = ub(W); xsq = FQ(W, hu(W))
    mid_ = fadd(FQ(D[2]), fmul(xq, FQ(D[3])))
    dp = long ? fadd(mid_, fmul(xsq, fadd(FQ(D[4]), fmul(xq, fadd(FQ(D[5]), fmul(xq, fadd(FQ(D[6]), fmul(xq, FQ(D[7]))))))))) : mid_
    tt1 = fmul(xq, dp)
    Xm = ub(X)
    @assert abs(D[1]) >= mag(tt1)
    t2m = Xm + Wm/2
    e_t2l = hu(hu(t2m) + big(2.0)^-54*Wm)
    t2l = hu(t2m) + big(2.0)^-54*Wm + e_t2l
    e_cube = mul122(Xm, Wm, hu(Wm))
    t1m = abs(D[1]) + mag(tt1)
    tt3m = Xm*Wm*t1m*(1 + big(2.0)^-50)
    @assert lb(X) - Wm/2 >= tt3m
    e = Xm*Wm*tt1.e + t1m*e_cube + tt3m*EM22 + e_t2l + add22(t2m, t2l, tt3m, hu(tt3m))
    deg = long ? 9 : 5
    P(s) = (q = Arb(D[deg - 2]) + 0*s; for k in deg-1:-1:3; q = Arb(D[k - 2]) + s*q; end; s + s*s/2 + s*s*s*q - (exp(s) - 1))
    A = tm(P, X, 8.0, 1.0; K = 14)
    (ub(A) + e)/lb(abs(expm1(X)))
end
# x scaled by 2^-s, then s steps expm1(2y) = (2 + expm1(y)) expm1(y)
function directX(a::BigFloat, b::BigFloat, s::Int, long::Bool)
    X = ival(a, b)
    EPS = rad(directpoly(X, long))
    for j in 0:s-1
        e = expm1(Arb(2)^j*X)
        r1l = hu(hu(ub(2 + e)) + hu(ub(abs(e))))/lb(2 + e)
        EPS = (1 + EPS)*(1 + e*EPS/(2 + e))*(1 + rad(r1l))*(1 + rad(EM22)) - 1
    end
    F = expm1(Arb(2)^s*X)
    metrics(lb(abs(F)), ub(abs(F)*EPS))
end

function run()
    @assert expmaj(big(1.04), Dict(k => D[k-2] for k in 3:9)) <= 8 && expmaj(big(1.04), Dict(k => D[k-2] for k in 3:5)) <= 8
    @printf("  expm1 common path: table 2^%.2f, |r| <= %.4e, reduction 2^%.2f, exp part relative error <= 2^%.3f\n",
            l2(ETAB), Float64(RMAX), l2(ERED), l2(EXPEPS))
    xdir = dbword(cint(H, "DIRECTINTERVALBOUND"), 0)          # direct path: |x| < xdir
    @assert xdir == 0.25
    cp = bnb(commonX, geopieces(big(0.25), big(OVF), 1 + big(2.0)^-6); relw = big(2.0)^-20, maxit = 600)
    cn = bnb(commonX, [(-b, -a) for (a, b) in geopieces(big(0.25), -big(MINUS1), 1 + big(2.0)^-6)]; relw = big(2.0)^-20, maxit = 600)
    xret = dbword(cint(H, "RETURNXBOUND"), 0)                 # |x| < xret returns x
    xlong = dbword(cint(H, "SPECIALINTERVALBOUND") + 1, 0)    # long polynomial from here on
    dw = [(big(0.0), big(0.0)) for _ in 1:2]
    for sg in (1, -1)
        for (lo, hi, long) in ((xret, xlong, false), (xlong, big(2.0)^-5, true))
            pcs = geopieces(lo, hi, 1 + big(2.0)^-4)
            sg < 0 && (pcs = [(-b, -a) for (a, b) in pcs])
            r = bnb((a, b) -> directX(a, b, 0, long), pcs; relw = big(2.0)^-16, maxit = 1500)
            dw = [supall(dw[o], r[o]) for o in 1:2]
        end
        for s in 1:3
            pcs = linpieces(big(2.0)^-6, big(2.0)^-5, 64)
            sg < 0 && (pcs = [(-b, -a) for (a, b) in pcs])
            r = bnb((a, b) -> directX(a, b, s, true), pcs; relw = big(2.0)^-16, maxit = 800)
            dw = [supall(dw[o], r[o]) for o in 1:2]
        end
    end
    [rowdir("expm1 common path x > 0 directed", cp[1][1], H, "ROUNDCSTCOMMONRD"),
     rowdir("expm1 common path x < 0 directed", cn[1][1], H, "ROUNDCSTCOMMONRD"),
     rowrn("expm1 common path RN relative", max(cp[2][1], cn[2][1]), H, "ROUNDCSTCOMMONRN"),
     rowdir("expm1 direct path directed", dw[1][1], H, "ROUNDCSTDIRECTRD"),
     rowrn("expm1 direct path RN relative", dw[2][1], H, "ROUNDCSTDIRECTRN")]
end
end
