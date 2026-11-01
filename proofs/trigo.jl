# trigo_fast.c case 2 (small arguments, no range reduction): sin, cos, tan

module TrigProof
include("common.jl")

const H = "trigo_fast.h"
const s3, s5, s7 = cdbnumber(H, "s3"), cdbnumber(H, "s5"), cdbnumber(H, "s7")
const c2, c4, c6 = cdbnumber(H, "c2"), cdbnumber(H, "c4"), cdbnumber(H, "c6")
const t3h, t3l, t5, t7, t9, t11 = (cdbnumber(H, n) for n in ("t3h", "t3l", "t5", "t7", "t9", "t11"))

# M bounds |p - f| on discs of radius rho around centres in [0, 1/16], from positive-coefficient majorants
majsin(R) = R + R^3*(abs(s3) + R^2*(abs(s5) + R^2*abs(s7))) + ub(sinh(Arb(R)))
majcos(R) = 1 + R^2*(abs(c2) + R^2*(abs(c4) + R^2*abs(c6))) + ub(cosh(Arb(R)))
majtan(R) = R + R^3*(abs(t3h) + abs(t3l) + R^2*(t5 + R^2*(t7 + R^2*(t9 + R^2*t11)))) + ub(tan(Arb(R)))
@assert majsin(big(1.0625)) <= 5 && majcos(big(1.0625)) <= 5 && majtan(big(0.5625)) <= 2

function sinpiece(X::Arb, rn::Bool)
    W = X*X; xx = FQ(W, hu(W))
    a = fmul(xx, FQ(s7)); a = fadd(FQ(s5), a); a = fmul(xx, a); a = fadd(FQ(s3), a)
    y = rn ? fmul(fmul(xx, a), FQ(X)) : fmul(fmul(FQ(X), xx), a)
    @assert lb(X) >= mag(y)*(1 + big(2.0)^-40)
    D = tm(s -> s + s*s*s*(Arb(s3) + s*s*(Arb(s5) + s*s*Arb(s7))) - sin(s), X, 5.0, 1.0)
    metrics(lb(sin(X)), ub(D) + y.e)
end
function cospiece(X::Arb)
    W = X*X; x2 = FQ(W, hu(W))
    a = fmul(x2, FQ(c6)); a = fadd(FQ(c4), a); a = fmul(x2, a); a = fadd(FQ(c2), a); tc = fmul(x2, a)
    @assert mag(tc) <= 1
    D = tm(s -> 1 + s*s*(Arb(c2) + s*s*(Arb(c4) + s*s*Arb(c6))) - cos(s), X, 5.0, 1.0)
    metrics(lb(cos(X)), ub(D) + tc.e)
end
function tanpiece(X::Arb)
    W = X*X; x2 = FQ(W, hu(W))
    a = fmul(x2, FQ(t11)); a = fadd(FQ(t9), a); a = fmul(x2, a); a = fadd(FQ(t7), a); a = fmul(x2, a); a = fadd(FQ(t5), a)
    b = fadd(FQ(t3l), fmul(x2, a)); tt = fmul(x2, fadd(FQ(t3h), b)); y = fmul(FQ(X), tt)
    @assert lb(X) >= mag(y)*(1 + big(2.0)^-40)
    D = tm(s -> s + s*s*s*(Arb(t3h) + Arb(t3l) + s*s*(Arb(t5) + s*s*(Arb(t7) + s*s*(Arb(t9) + s*s*Arb(t11))))) - tan(s), X, 2.0, 0.5)
    metrics(lb(tan(X)), ub(D) + y.e)
end
# rndcst of tan_rn for |x| with high word hi
function tan_rndcst(x::Float64)
    hi = Int64(reinterpret(UInt64, abs(x)) >> 32)
    dbword(0x3ff00000 + (((hi & 0x000fffff) + 0x00100000) >> (0x3ff + 2 - (hi >> 20))), 0xffffffff)
end

function run()
    lo_s = dbword(cint(H, "XMAX_RETURN_X_FOR_SIN"), 0); hi_s = prevfloat(dbword(cint(H, "XMAX_SIN_CASE2"), 0))
    lo_cn = dbword(cint(H, "XMAX_RETURN_1_FOR_COS_RN"), 0); lo_cd = dbword(cint(H, "XMAX_RETURN_1_FOR_COS_RDIR"), 0)
    hi_c = prevfloat(dbword(cint(H, "XMAX_COS_CASE2"), 0))
    lo_t = dbword(cint(H, "XMAX_RETURN_X_FOR_TAN"), 0); hi_t = prevfloat(dbword(cint(H, "XMAX_TAN_CASE2"), 0))
    @assert hi_t < 1/16 && hi_s < 1/16
    B(f, lo, hi) = bnb((a, b) -> f(ival(a, b)), geopieces(lo, hi, 1 + big(2.0)^-5); relw = big(2.0)^-16, maxit = 3000)
    sd = B(X -> sinpiece(X, false), lo_s, hi_s); sn = B(X -> sinpiece(X, true), lo_s, hi_s)
    cd = B(cospiece, lo_cd, hi_c); cn = B(cospiece, lo_cn, hi_c)
    td = B(tanpiece, lo_t, hi_t)
    # tan_rn: relative bound against the eps allowed by the constant at the bottom of each piece
    worst = big(0.0)
    for (a, b) in geopieces(lo_t, hi_t, 1 + big(2.0)^-8)
        q = tanpiece(ival(a, b))[2]/eps_rn(tan_rndcst(Float64(a, RoundDown)))
        if q > 1
            for (aa, bb) in linpieces(a, b, 64)
                worst = max(worst, tanpiece(ival(aa, bb))[2]/eps_rn(tan_rndcst(Float64(aa, RoundDown))))
            end
        else
            worst = max(worst, q)
        end
    end
    [rowdir("sin case 2 directed", sd[1][1], H, "EPS_SIN_CASE2"),
     rowrn("sin case 2 RN relative", sn[2][1], H, "RN_CST_SIN_CASE2"),
     rowdir("cos case 2 directed", cd[1][1], H, "EPS_COS_CASE2"),
     rowrn("cos case 2 RN relative", cn[2][1], H, "RN_CST_COS_CASE2"),
     rowdir("tan case 2 directed", td[1][1], H, "EPS_TAN_CASE2"),
     Row("tan_rn case 2 relative / eps of rndcst", worst, big(1.0), "$(ctext("trigo_fast.c", "rndcst.i[HI] = 0x3ff00000")) rndcst")]
end
end
