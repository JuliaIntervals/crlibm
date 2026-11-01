# atan_fast.c quick phase (atan_quick), used by atan and atanpi

module AtanProof
include("common.jl")

const H = "atan_fast.h"
const TA = reshape(cdbnumbers(H, "arctan_table"), 4, 62)      # a[i], b[i], atan_b[i] hi, lo
const EPS = cdoubles(H, "epsilon"); const RNC = cdoubles(H, "rncst")
const c9, c7, c5, c3 = cdoubles(H, "coef_poly")
const MINRED = cdouble(H, "MIN_REDUCTION_NEEDED")
const INVPIH, INVPIL = cdouble(H, "INVPIH"), cdouble(H, "INVPIL")
@assert abs(INVPIL) <= hu(INVPIH) && issorted(TA[1, :]) && MINRED >= TA[1, 1]
const PI = Arb(π)
const ED = big(2.0)^-102     # Div22 relative error
const EADD = big(2.0)^-104   # Add22 relative error in 1 + x*b[i]
a_(i) = TA[1, i+1]; b_(i) = TA[2, i+1]; bhi(i) = TA[3, i+1]; blo(i) = TA[4, i+1]

# index i chosen by atan_quick for x > MIN_REDUCTION_NEEDED
function findi(x)
    x > b_(61) && return 61
    i = 31
    for s in (16, 8, 4, 2, 1); i = x < a_(i) ? i - s : i + s; end
    x < a_(i) && (i -= 1)
    i
end

const NA = 60
const DC = [A0; [Arb(k == 1 ? c3 : k == 2 ? c5 : k == 3 ? c7 : k == 4 ? c9 : 0.0) - Arb((-1)^k)/(2k + 1) for k in 1:NA]]
Dser(w) = (s = 0*w + DC[NA+1]; for k in NA-1:-1:1; s = s*w + DC[k+1]; end; s*w)
# |t*Dser(t^2)| <= 2 on discs of radius 1/2 around centres |t| <= 0.0127
@assert (R = big(0.5127); R*sum(abs(DC[k+1])*R^(2k) for k in 1:NA) <= 2)
# enclosure of t + t*Q(t^2) - atan(t) over T
function approx(T::Arb)
    @assert ub(mid(T)) <= 0.0127
    tmax = ub(T)
    tm(s -> s*Dser(s*s), T, 2.0, 0.5; K = 12) + rad(tmax^(2NA + 3)/(1 - tmax^2))
end
# rounding error of Xredhi*q and its magnitude
function polyround(T::Arb)
    W = T*T; x2 = FQ(W, hu(W))
    a = fmul(x2, FQ(c9)); a = fadd(FQ(c7), a); a = fmul(x2, a); a = fadd(FQ(c5), a); a = fmul(x2, a); a = fadd(FQ(c3), a)
    m = fmul(FQ(T), fmul(x2, a))
    m.e, mag(m)
end
# atan and atanpi: (directed metric / epsilon, relative / eps allowed by rncst)
function ratios(F::Arb, e, ide)
    rho = ub(Arb(INVPIH) + Arb(INVPIL) - 1/PI); Fm = ub(F)
    ep = (abs(INVPIH) + abs(INVPIL))*e + Fm*rho + Fm*0.3184*EM22   # Mul22(INVPIH, INVPIL, atanhi, atanlo)
    m1 = metrics(lb(F), e); m2 = metrics(lb(F/PI), ep)
    [m1[1]/EPS[ide+1], m1[2]/eps_rn(RNC[ide+1]), m2[1]/EPS[ide+1], m2[2]/eps_rn(RNC[ide+1])]
end

function reduced(i, a::BigFloat, b::BigFloat)
    bi = b_(i)
    X = ival(a, b)
    # x0 = 1 + x*b[i] by Mul12 and Add22 (exponent of first operand >= the other's), then Div22
    P = X*Arb(bi); pm = ub(P)*(1 + big(2.0)^-52)
    if a > 1
        @assert lb(P) >= 1
        d0 = add22(pm, hu(pm), big(1.0), big(0.0))
    else
        @assert pm < 2
        d0 = add22(big(1.0), big(0.0), pm, hu(pm))
    end
    x0 = 1 + P; d0 /= lb(x0)
    num = X - Arb(bi)
    y0lo = lb(x0)*(1 - big(2.0)^-52); y0m = ub(x0)*(1 + big(2.0)^-52)
    dD = div22rel(y0lo, y0m, i == 61)
    @assert (dD + d0)/(1 - d0) <= ED + EADD
    XR = num/(1 + X*Arb(bi))
    xm = ub(XR)
    T = XR*(1 + rad(ED + EADD)); T = T + rad(hu(ub(T)))      # ball of Xredhi
    tmax = ub(T); hlo = hu(tmax)                             # |Xredlo|
    eA = ub(approx(T))
    eLo = hlo*tmax^2*(1 + big(2.0)^-40) + hlo^2              # Xredlo against atan(Xredhi + Xredlo) - atan(Xredhi)
    eX = xm*(ED + EADD)*(1 + big(2.0)^-40)                   # atan of the computed Xred against atan(Xred)
    tau = ub(Arb(bhi(i)) + Arb(blo(i)) - atan(Arb(bi)))
    eq, mq = polyround(T)
    l1 = hlo + abs(blo(i)); e1 = hu(l1)
    l2v = l1 + e1 + mq + eq; e2 = hu(l2v)
    @assert abs(bhi(i)) >= tmax
    tl = hu(abs(bhi(i)) + tmax)
    s = tl + l2v + e2; e3 = hu(s)
    e = eA + eLo + eX + tau + eq + e1 + e2 + e3
    ratios(atan(X), e, i < 10 ? 0 : 1)
end
function noreduction(a::BigFloat, b::BigFloat)
    X = ival(a, b)
    eq, mq = polyround(X)
    @assert lb(X) >= mq
    ratios(atan(X), ub(approx(X)) + eq, 2)
end

# relative error of Div22 for a numerator xh (+ xl, |xl| <= 2^-53|xh|, if withxl) and a
# normalized denominator yh + yl with yhlo <= |yh| <= yhm; hu(v) <= 2^-53 v makes it scale-free
function div22rel(yhlo, yhm, withxl)
    u = big(2.0)^-53
    yl = u*yhm
    m3 = yhm*u*(withxl ? 2(1 + 2u) : 1)
    mp = (1 + 2u)*yl
    e3 = u*m3; ep = u*mp
    m4 = m3 + e3 + mp + ep; e4 = u*m4
    ((m3 + mp)*yl/(yhlo*(yhlo - yl)) + (e3 + ep + e4)/yhlo + u*m4/yhlo)/(1 - 4u)
end

# x - b[i] is exact on [lo, hi]: both are multiples of 2^(emin-52) and |x - b| < 2^(emin+1)
function exactsub(lo, hi, bi)
    for ex in exponent(lo):exponent(hi)
        xl = max(big(lo), big(2.0)^ex); xh = min(big(hi), big(2.0)^(ex+1))
        emin = min(ex, exponent(bi))
        max(abs(xl - bi), abs(xh - bi)) < big(2.0)^(emin+1) || return false
    end
    true
end

function run()
    pts = sort(unique([nextfloat(MINRED); TA[1, :]; prevfloat.(TA[1, :]); b_(61); nextfloat(b_(61)); prevfloat(2.0^54)]))
    ranges = Dict{Int,Vector{Float64}}()
    for p in pts
        p <= MINRED && continue
        push!(get!(ranges, findi(p), Float64[]), p)
    end
    @assert sort(collect(keys(ranges))) == 0:61
    names = ["atan directed / epsilon", "atan RN relative / eps of rncst", "atanpi directed / epsilon", "atanpi RN relative / eps of rncst"]
    w = Dict(ide => [(big(0.0), big(0.0)) for _ in 1:4] for ide in 0:2)
    for i in 0:61
        lo = minimum(ranges[i]); hi = i == 61 ? prevfloat(2.0^54) : maximum(ranges[i])
        @assert i == 61 || findi(nextfloat(hi)) == i + 1
        @assert exactsub(lo, i == 61 ? b_(61) : hi, b_(i))  # Add12 is used above b[61]
        pcs = Tuple{BigFloat,BigFloat}[]
        segs = lo < 1 < hi ? [(lo, 1.0), (nextfloat(1.0), hi)] : [(lo, hi)]
        for (l, h) in segs
            append!(pcs, i == 61 ? geopieces(l, h, 1 + big(2.0)^-6) : linpieces(l, h, 128 ÷ length(segs)))
        end
        res = bnb((a, b) -> reduced(i, a, b), pcs; relw = big(2.0)^-18, maxit = 800)
        ide = i < 10 ? 0 : 1
        w[ide] = [supall(w[ide][o], res[o]) for o in 1:4]
    end
    lo = dbword(0x3E400000, 0)
    ctext("atan_fast.c", "absxhi < 0x3E400000")
    res = bnb(noreduction, geopieces(lo, MINRED, 1 + big(2.0)^-6); relw = big(2.0)^-18, maxit = 2000)
    w[2] = [supall(w[2][o], res[o]) for o in 1:4]
    rows = Row[]
    for (fn, o) in (("atan", 0), ("atanpi", 2)), ide in 0:2
        push!(rows, Row("$fn directed, $(("i < 10", "i >= 10", "no reduction")[ide+1])", w[ide][o+1][1]*EPS[ide+1], EPS[ide+1],
                        "$(ctext(H, "static const double epsilon")) epsilon[$ide]"))
    end
    for (fn, o) in (("atan", 1), ("atanpi", 3)), ide in 0:2
        e = eps_rn(RNC[ide+1])
        push!(rows, Row("$fn RN relative, $(("i < 10", "i >= 10", "no reduction")[ide+1])", w[ide][o+1][1]*e, e,
                        "$(ctext(H, "static const double rncst")) rncst[$ide] = $(RNC[ide+1]), allows"))
    end
    rows
end
end
