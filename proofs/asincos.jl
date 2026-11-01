# asincos.c quick phases: asin, acos, asinpi, acospi (paths 1, 2 and 3),
# and sqrt12_64_unfiltered (crlibm_private.h) used by path 3

module AsinProof
include("common.jl")

const H = "asincos.h"
const P = "crlibm_private.h"
const SC = [cdouble(P, "SQRTPOLYC$k") for k in 0:4]
const P9 = [cdouble(H, "p9_quick_coeff_$(k)h") for k in 0:11]
const P0 = Dict(k => cdouble(H, "p0_quick_coeff_$(k)h") for k in 3:2:19)
const TQ = cdoubles(H, "asinPolyQuickMiddleTbl")
const MI9 = cdouble(H, "MI_9")
const PIHALFH, PIHALFM = cdouble(H, "PIHALFH"), cdouble(H, "PIHALFM")
const PIH, PIM = cdouble(H, "PIH"), cdouble(H, "PIM")
const RPH, RPM = cdouble(H, "RECPRPIH"), cdouble(H, "RECPRPIM")
@assert cdouble(H, "MRECPRPIH") == -RPH && cdouble(H, "MRECPRPIM") == -RPM
@assert all(abs(l) <= hu(abs(h)) for (h, l) in ((PIHALFH, PIHALFM), (PIH, PIM), (RPH, RPM)))
@assert length(TQ) == 128
const PI = Arb(π)
const RHO = Arb(RPH) + Arb(RPM) - 1/PI
const U = big(2.0)^-53

# ---------------------------------------------------------------- input ranges of the paths

# index in asin_*: bits 16..19 of the high word of 1.0 + x*x
pathindex(x::Float64) = Int((reinterpret(UInt64, 1.0 + x*x) >> 48) & 0xf)
# smallest positive double x with pathindex(x) >= t (pathindex is monotone in x)
function firstx(t)
    lo, hi = reinterpret(UInt64, 0.1), reinterpret(UInt64, 1.0)
    while hi - lo > 1
        m = (lo + hi) >> 1
        pathindex(reinterpret(Float64, m)) >= t ? (hi = m) : (lo = m)
    end
    reinterpret(Float64, hi)
end
const XB = [firstx(t) for t in 1:9]     # path 2 interval j: [XB[j+1], XB[j+2]); path 3: [XB[9], 1)
const XA2 = dbword(cint(H, "EXTRABOUND2"), 0)       # p0_quick branch a below this
const XA = dbword(cint(H, "EXTRABOUND") + 1, 0)     # branch c from this on
const XLOW = Dict(:asin => dbword(cint(H, "ASINSIMPLEBOUND"), 0), :acos => dbword(cint(H, "ACOSSIMPLEBOUND"), 0),
                  :asinpi => dbword(cint(H, "ASINPISIMPLEBOUND"), 0), :acospi => dbword(cint(H, "ACOSPISIMPLEBOUND"), 0))

# ---------------------------------------------------------------- sqrt12_64_unfiltered

r0s(m) = Arb(SC[1]) + m*(Arb(SC[2]) + m*(Arb(SC[3]) + m*(Arb(SC[4]) + m*Arb(SC[5]))))
function r0FQ(M::Arb)
    a = fmul(FQ(M), FQ(SC[5])); a = fadd(FQ(SC[4]), a)
    a = fmul(FQ(M), a); a = fadd(FQ(SC[3]), a)
    a = fmul(FQ(M), a); a = fadd(FQ(SC[2]), a)
    a = fmul(FQ(M), a); fadd(FQ(SC[1]), a)
end
# Newton step r' = r(3 - m r^2)/2: if r = (1+d)/sqrt(m) then r' = (1+F1(d))/sqrt(m)
F1(d) = -Arb(3)/2*d*d - d*d*d/2
# relative rounding error of r' = (0.5*r)*(3 - m*(r*r)) in double, when |d| <= dm
eta(dm) = (q = (1 + dm)^2; (1 + U)^2*(1 + q*(2U + U^2)/(3 - q)) - 1)
function delta0(M::Arb)
    s = ArbSeries((M, 1), degree = 1)
    D = (r0s(s)*sqrt(s) - 1)[1]
    c = mid(M)
    r0s(c)*sqrt(c) - 1 + D*rad(radius(M)) + rad(r0FQ(M).e*ub(sqrt(M)))
end
# rounding of the double-double step and of Mul122(m, r3), relative to sqrt(m)
function ddstep(M::Arb, D2::Arb)
    u = big(2.0)^-52
    mm = ub(M); R2 = (1 + D2)/sqrt(M); r2 = ub(R2)
    sq = r2^2*(1 + u); mr = mm*r2*(1 + u); ph = 1.5*r2*(1 + u)
    Pm = mr*sq*(1 + u)
    eM = mul22(mr, hu(mr), sq, hu(sq))
    @assert 1.5*lb(R2) >= Pm/2*(1 + u)
    eA = add22(ph, hu(ph), Pm/2, hu(Pm/2) + eM/2)
    r3 = ph + Pm/2
    (mm*(eM/2 + eA) + mul122(mm, r3, hu(r3)))/lb(sqrt(M))
end
function deltas(M::Arb)
    d0 = delta0(M)
    d1 = (1 + F1(d0))*(1 + rad(eta(ub(d0)))) - 1
    d2 = (1 + F1(d1))*(1 + rad(eta(ub(d1)))) - 1
    @assert ddstep(M, d2) <= big(2.0)^-99
    F1(d2) + rad(big(2.0)^-99)
end
# enclosures of the relative error of sqrt12_64_unfiltered on 2^14 pieces of m in [1/2, 2]
const NM = 1 << 14
const SQ = let t = Vector{Tuple{BigFloat,BigFloat}}(undef, NM)
    for j in 1:NM
        d = deltas(ival(big(0.5) + big(1.5)*(j - 1)/NM, big(0.5) + big(1.5)*j/NM))
        t[j] = (BigFloat(Arblib.lbound(Arf, d)), BigFloat(Arblib.ubound(Arf, d)))
    end
    t
end
function dsq(mlo, mhi)
    j1 = clamp(floor(Int, (mlo - 0.5)/1.5*NM) + 1, 1, NM); j2 = clamp(floor(Int, (mhi - 0.5)/1.5*NM) + 1, 1, NM)
    Arb((minimum(SQ[j][1] for j in j1:j2), maximum(SQ[j][2] for j in j1:j2)))
end

# ---------------------------------------------------------------- path 3: p9_quick, |x| >= 0.75

p9ex(z) = (s = Arb(P9[12]) + 0*z; for k in 11:-1:1; s = Arb(P9[k]) + z*s; end; s)
# g = acos(x)/sqrt(2(1-x)) = asin(s)/s with s^2 = w = (1-x)/2 = sum_n HC[n] w^n, HC[n] <= 1
const HC = [Arb(binomial(big(2n), big(n)))/(Arb(4)^n*(2n + 1)) for n in 0:200]
hser(w, N) = (s = HC[N+1] + 0*w; for n in N-1:-1:0; s = s*w + HC[n+1]; end; s)
wof(z) = (1 - Arb(MI9) - z)/2
gfun(Z::Arb) = (wm = ub(wof(Z)); @assert wm < 0.2; hser(wof(Z), 90) + rad(wm^91/(1 - wm)))
# |p9 - g_90| <= 2.67 on discs of radius 1 around centres |z| <= 0.126
@assert (R = big(1.126); Wd = (abs(1 - big(MI9)) + R)/2;
         sum(abs(P9[k+1])*R^k for k in 0:11) + ub(asin(sqrt(Arb(Wd)))/sqrt(Arb(Wd)))) <= 2.67
function eapp(Z::Arb)
    r = radius(Z); @assert r < 0.01 && ub(mid(Z)) <= 0.126
    E = tm(s -> p9ex(s) - hser(wof(s), 90), Z, 2.67, 1.0; K = 14)
    wm = ub(wof(Z))
    (E + rad(wm^91/(1 - wm)))/gfun(Z)
end
function p9round(Z::Arb)
    zq = FQ(Z)
    a = FQ(P9[12])
    for k in 11:-1:3; a = fadd(FQ(P9[k]), fmul(a, zq)); end
    p20 = fmul(a, zq)
    @assert abs(P9[2]) >= mag(p20)
    t21 = Arb(P9[2]) + p20.v; t21m = hu(ub(t21) + p20.e)
    zm = ub(Z); t21h = ub(t21) + p20.e
    t3 = zm*t21m; t2 = hu(zm*t21h)
    eM = hu(t3) + hu(t2 + t3 + hu(t3))
    t22 = zm*t21h*(1 + big(2.0)^-52) + t2 + t3
    @assert abs(P9[1]) >= t22
    eA = hu(hu(abs(P9[1]) + t22) + hu(t22) + t2 + t3)
    zm*p20.e + eM + eA
end
# zp = 2(1-|x|) in [a, b] inside the binade [2^k, 2^(k+1)]
function piece3(a::BigFloat, b::BigFloat, k::Int)
    ZP = ival(a, b)
    Z = 1 - ZP/2 - Arb(MI9)                     # z = |x| - MI_9 and zp are exact (Sterbenz)
    mlo, mhi = isodd(k) ? (a/big(2.0)^(k+1), b/big(2.0)^(k+1)) : (a/big(2.0)^k, b/big(2.0)^k)
    ep = eapp(Z) + rad(p9round(Z)/lb(gfun(Z)))
    eps1 = (1 + dsq(mlo, mhi))*(1 + ep)*(1 + rad(EM22)) - 1
    AC = acos(1 - ZP/2)
    eps2 = (1 + eps1)*(1 + RHO*PI)*(1 + rad(EM22)) - 1
    T1m = ub(AC)*(1 + big(2.0)^-50); T2m = T1m/3
    e1 = ub(AC*eps1); e2 = ub(AC/PI*eps2)
    [metrics(lb(AC), e1);                                                     # acos, x > 0
     metrics(lb(PI - AC), e1 + add22(PIH, hu(PIH), T1m, hu(T1m)) + ub(Arb(PIH) + Arb(PIM) - PI));        # acos, x < 0
     metrics(lb(PI/2 - AC), e1 + add22(PIHALFH, hu(PIHALFH), T1m, hu(T1m)) + ub(Arb(PIHALFH) + Arb(PIHALFM) - PI/2)); # asin
     metrics(lb(AC/PI), e2);                                                  # acospi, x > 0
     metrics(lb(1 - AC/PI), e2 + add122(1, T2m, hu(T2m)));                    # acospi, x < 0
     metrics(lb(0.5 - AC/PI), e2 + add122(0.5, T2m, hu(T2m)))]                # asinpi
end

# ---------------------------------------------------------------- paths 1 and 2

# metrics of the six results, given asin(|x|) in AS and the relative error EPS of the computed asin
function compose(AS::Arb, EPS::Arb, path)
    asm = ub(AS)*(1 + big(2.0)^-50)
    e = ub(AS*EPS)
    eps2 = (1 + EPS)*(1 + RHO*PI)*(1 + rad(EM22)) - 1
    eA = (path == 1 ? add22 : add22cond)(PIHALFH, hu(PIHALFH), asm, hu(asm)) + ub(Arb(PIHALFH) + Arb(PIHALFM) - PI/2)
    e2 = ub(AS/PI*eps2) + add122(0.5, asm/3, hu(asm/3))
    [metrics(lb(AS), e);                     # asin
     metrics(lb(AS/PI), ub(AS/PI*eps2));     # asinpi
     metrics(lb(PI/2 - AS), e + eA);         # acos, x > 0
     metrics(lb(PI/2 + AS), e + eA);         # acos, x < 0
     metrics(lb(0.5 - AS/PI), e2);           # acospi, x > 0
     metrics(lb(0.5 + AS/PI), e2)]           # acospi, x < 0
end
const NAMES = ["asin", "asinpi", "acos+", "acos-", "acospi+", "acospi-"]

# |q(w) - h_N(w)| <= 3.5 on discs of radius 1/2 around centres w <= 1/16
@assert (R = big(0.5625); 1 + sum(abs(P0[k])*R^((k-1)÷2) for k in 3:2:19) + ub(asin(sqrt(Arb(R)))/sqrt(Arb(R)))) <= 3.5
function p0piece(X::Arb, branch)
    W = X*X; Wm = ub(W); Xm = ub(X); Xl = lb(X)
    if branch == :a
        t2 = fmul(FQ(P0[3]), FQ(X)); t1 = fmul(FQ(X), FQ(X)); t3 = fmul(t1, t2)
        @assert Xl >= mag(t3)
        Er = t3.e
        q = [1.0, P0[3]]
    else
        wq = FQ(W, hu(W))                            # x2h = RN(x^2), x2h + x2m = x^2
        if branch == :c
            a = FQ(P0[19])
            for k in 17:-2:7; a = fadd(FQ(P0[k]), fmul(a, wq)); end
            a = fadd(FQ(P0[5]), fmul(a, wq))
            q = [1.0; [P0[k] for k in 3:2:19]]
        else
            a = FQ(P0[5]); q = [1.0, P0[3], P0[5]]
        end
        p16 = fmul(a, wq)
        @assert abs(P0[3]) >= mag(p16)
        p17 = Arb(P0[3]) + p16.v
        t3 = Xm*hu(Wm); t2 = hu(Xm*Wm)
        e18 = hu(t3) + hu(t2 + t3 + hu(t3))           # Mul122(x, x2h, x2m)
        X3 = Xm*Wm
        p19m = (ub(p17) + p16.e)*(X3 + e18)*(1 + EM22)
        @assert Xl >= p19m*(1 + big(2.0)^-50)
        Er = X3*p16.e + (ub(p17) + p16.e)*e18 + p19m*EM22 + add122(Xm, p19m, hu(p19m))
    end
    N = 120; K = 16
    wmid = mid(W); r = radius(W); @assert r < 0.05 && ub(wmid) <= 1/16
    sw = ArbSeries((wmid, 1), degree = K)
    Ds = ArbSeries((A0,), degree = K)
    for n in N:-1:1
        Ds = Ds*sw + ((n < length(q) ? Arb(q[n+1]) : A0) - HC[n+1])
    end
    Ds = Ds*sw
    t = rad(r); D = A0
    for k in K:-1:0; D = D*t + Ds[k]; end
    D += rad(big(3.5)*(r/big(0.5))^(K + 1)/(1 - r/big(0.5))) + rad(Wm^(N + 1)/(1 - Wm))
    Hs = hser(W, N) + rad(Wm^(N + 1)/(1 - Wm))
    compose(X*Hs, D/Hs + rad(Er/(Xl*lb(Hs))), 1)
end

function p2piece(j::Int, X::Arb)
    c = TQ[16j+1:16j+16]; mi = c[1]
    c0h, c0m, c1h, c1m = c[2], c[3], c[4], c[5]
    ck = Dict(k => c[k+4] for k in 2:12)
    @assert lb(X) >= mi/2 && ub(X) <= 2mi            # z = |x| - mi exact (Sterbenz)
    Z = X - Arb(mi); Zm = ub(Z)
    a = FQ(ck[12])
    for k in 11:-1:3; a = fadd(FQ(ck[k]), fmul(a, FQ(Z))); end
    p20 = fmul(a, FQ(Z))
    @assert abs(ck[2]) >= mag(p20)
    t21 = abs(ck[2]) + mag(p20); E21 = p20.e
    function madd(ch, cl, bmag, Eb)                   # MulAdd212(ch, cl, z, bh, bl), |b| <= bmag
        t1 = Zm*bmag; @assert (abs(ch) - t1) >= t1*(1 + big(2.0)^-50)
        t2 = hu(t1); t3 = abs(ch) + t1; t4 = hu(t3)
        t5 = Zm*hu(bmag); e5 = hu(t5)
        t6 = abs(cl) + t2; e6 = hu(t6)
        t7 = t5 + e5 + t6 + e6; e7 = hu(t7)
        t8 = t7 + e7 + t4; e8 = hu(t8)
        @assert (abs(ch) - t1)*(1 - big(2.0)^-50) >= t8
        (t3 + t8, Zm*Eb + e5 + e6 + e7 + e8)
    end
    t22, E22 = madd(c1h, c1m, t21, E21)
    t23, E23 = madd(c0h, c0m, t22, E22)
    r = radius(X); cx = ub(mid(X)); @assert r < 2.0^-10 && cx + 0.2 < 1
    cs = [big(c0h) + c0m; big(c1h) + c1m; [big(ck[k]) for k in 2:12]]
    Mc = sum(abs(cs[k+1])*(0.2 + abs(cx - mi))^k for k in 0:12) + ub(asin(Arb(cx + 0.2)))
    Pz(s) = (zs = s - Arb(mi); v = Arb(ck[12]) + 0*zs;
             for k in 11:-1:2; v = Arb(ck[k]) + zs*v; end;
             Arb(c0h) + Arb(c0m) + zs*(Arb(c1h) + Arb(c1m) + zs*v))
    v = tm(s -> Pz(s) - asin(s), X, Mc, 0.2)
    AS = asin(X)
    compose(AS, v/AS + rad(E23/lb(AS)), 2)
end

maskfor(names) = [n in names for n in NAMES for _ in 1:2]

function run()
    # sqrt12_64_unfiltered
    slo = minimum(t[1] for t in SQ); shi = maximum(t[2] for t in SQ)
    @printf("  sqrt12_64_unfiltered: relative error in [-2^%.3f, +2^%.2f]\n", l2(-slo), l2(shi))
    rows = [Row("sqrt12_64_unfiltered |relative error|", max(-slo, shi), big(2.0)^-62.4, "$(ctext(P, "2^(-62.4)")) comment 2^-62.4")]

    # path 3
    mx = big(0.0); mr = big(0.0)
    zlo = XB[9] - big(MI9); zhi = 1 - big(2.0)^-53 - MI9
    for (a, b) in linpieces(zlo, zhi, 4096)
        Z = ival(a, b); mx = max(mx, ub(eapp(Z))); mr = max(mr, p9round(Z)/lb(gfun(Z)))
    end
    @printf("  asin path 3: p9_quick approximation <= 2^%.3f, rounding <= 2^%.3f (relative)\n", l2(mx), l2(mr))
    zpmax = 2*(1 - big(XB[9]))
    p3 = Tuple{BigFloat,BigFloat,Int}[]
    w3 = [(big(0.0), big(0.0)) for _ in 1:12]
    for k in -53:-1
        lo = big(2.0)^k; hi = k == -1 ? big(0.5)*(1 + big(2.0)^-40) : big(2.0)^(k+1)
        pcs = linpieces(lo, hi, 32)
        append!(p3, [(a, b, k) for (a, b) in pcs])
        res = bnb((a, b) -> piece3(a, b, k), pcs; relw = big(2.0)^-14, maxit = 600)
        w3 = [supall(w3[o], res[o]) for o in 1:12]
    end
    @assert covers([(a, b) for (a, b, _) in p3], big(2.0)^-52, zpmax)
    d3(o) = w3[o][1]
    # acos+ acos- asin acospi+ acospi- asinpi, (directed, relative)
    append!(rows, [
        rowdir("asin   path 3 directed", d3(5), H, "RDROUNDCST9"),
        rowdir("acos   path 3 directed", max(d3(1), d3(3)), H, "RDROUNDCST9"),
        rowdir("asinpi path 3 directed", d3(11), H, "RDROUNDCST9"),
        rowdir("acospi path 3 directed", max(d3(7), d3(9)), H, "RDROUNDCST9"),
        rowrn("asin   path 3 RN relative", d3(6), H, "RNROUNDCST9"),
        rowrn("acos   path 3 RN relative", max(d3(2), d3(4)), H, "RNROUNDCST9"),
        rowrn("asinpi path 3 RN relative", d3(12), H, "RNROUNDCST9"),
        rowrn("acospi path 3 RN relative", max(d3(8), d3(10)), H, "RNROUNDCST9")])

    # paths 1 and 2
    w = [(big(0.0), big(0.0)) for _ in 1:12]
    xlow = minimum(values(XLOW))
    parts = [(:a, xlow, XLOW[:asin], maskfor(("asinpi", "acospi+", "acospi-"))),
             (:a, XLOW[:asin], XA2, trues(12)), (:b, XA2, XA, trues(12)), (:c, big(2.0)^-10, big(0.25)*(1 + big(2.0)^-50), trues(12))]
    @assert XLOW[:asin] <= XLOW[:acos] && XLOW[:acospi] <= XLOW[:asin] && XA <= XB[1]
    @assert big(prevfloat(XB[1])) <= big(0.25)*(1 + big(2.0)^-50)
    for (br, lo, hi, mask) in parts
        pcs = geopieces(lo, hi, 1 + big(2.0)^-4)
        res = bnb((a, b) -> [m ? v : big(0.0) for (v, m) in zip(p0piece(ival(a, b), br), mask)], pcs; relw = big(2.0)^-16, maxit = 1500)
        w = [supall(w[o], res[o]) for o in 1:12]
    end
    for j in 0:7
        lo = sqrt(big(j + 1)/16)*(1 - big(2.0)^-50); hi = sqrt(big(j + 2)/16)*(1 + big(2.0)^-50)
        @assert lo <= XB[j+1] && hi >= prevfloat(XB[j+2])
        res = bnb((a, b) -> p2piece(j, ival(a, b)), linpieces(lo, hi, 64); relw = big(2.0)^-16, maxit = 1500)
        w = [supall(w[o], res[o]) for o in 1:12]
    end
    d(o) = w[o][1]
    append!(rows, [
        rowdir("asin   paths 1-2 directed", d(1), H, "RDROUNDCST"),
        rowdir("acos   paths 1-2 directed", max(d(5), d(7)), H, "RDROUNDCST"),
        rowdir("asinpi paths 1-2 directed", d(3), H, "RDROUNDCST"),
        rowdir("acospi paths 1-2 directed", max(d(9), d(11)), H, "RDROUNDCST"),
        rowrn("asin   paths 1-2 RN relative", d(2), H, "RNROUNDCST"),
        rowrn("acos   paths 1-2 RN relative", max(d(6), d(8)), H, "RNROUNDCST"),
        rowrn("asinpi paths 1-2 RN relative", d(4), H, "RNROUNDCST"),
        rowrn("acospi paths 1-2 RN relative", max(d(10), d(12)), H, "RNROUNDCST")])
    rows
end
end
