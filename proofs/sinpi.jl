# trigpi.c: sinpiquick and the first-step test of sinpi_rn (rh == rh + 1.00001*rm)

module SinpiProof
include("common.jl")

const F = "trigpi.c"
const c1h, c1m, c3, c5, c7 = (cdouble(F, "sinpiquick_coeff_$(k)") for k in ("1h", "1m", "3h", "5h", "7h"))
const k2, k4, k6 = (cdouble(F, "cospiquick_coeff_$(k)h") for k in 2:2:6)
@assert cdouble(F, "cospiquick_coeff_0h") == 1.0 && abs(c1m) <= hu(c1h)
const ST = reshape(cdoubles("trigpi.h", "sincosTable"), 6, 64)    # sh, ch, sm, cm, sl, cl
@assert all(abs(ST[3, i]) <= hu(abs(ST[1, i])) && abs(ST[4, i]) <= hu(abs(ST[2, i])) for i in 1:64)
const RNC = parse(Float64, cfind(F, r"rh\s*==\s*rh\s*\+\s*([0-9.]+)\s*\*\s*rm")[1])
const ERN = eps_rn(RNC)
const EMA = big(2.0)^-100; const EM122 = big(2.0)^-102
const PI = Arb(π)
const R = big(1.0) + 1/big(256)
@assert R*(abs(c1h) + abs(c1m) + R^2*(abs(c3) + R^2*(abs(c5) + R^2*abs(c7)))) + ub(sinh(PI*R)) <= 30
@assert 1 + R^2*(abs(k2) + R^2*(abs(k4) + R^2*abs(k6))) + ub(cosh(PI*R)) <= 30

# errors of sy = sin(pi y), cy = cos(pi y), and bounds of |sy|, |cy|
function scpiece(Y::Arb)
    W = Y*Y; Wm = ub(W); Ym = ub(Y)
    x2h = FQ(W, hu(W))
    t4 = fmul(fadd(FQ(c5), fmul(FQ(c7), x2h)), x2h)
    t5 = abs(c3) + mag(t4)
    @assert abs(c3) >= mag(t4)
    t6m = abs(c1h) + Wm*t5
    @assert muladd22(abs(c1h), abs(c1m), Wm*(1 + big(2.0)^-52), hu(Wm), t5*(1 + big(2.0)^-52), hu(t5)) <= EMA*t6m
    @assert mul122(Ym, t6m*(1 + big(2.0)^-50), hu(t6m*(1 + big(2.0)^-50))) <= EM122*Ym*t6m
    esy = Ym*(Wm*t4.e + t6m*EMA) + Ym*t6m*EM122*(1 + big(2.0)^-50)
    u4 = fmul(fadd(FQ(k4), fmul(FQ(k6), x2h)), x2h)
    @assert abs(k2) >= mag(u4)
    t6c = Wm*(abs(k2) + mag(u4))
    ecy = Wm*u4.e + t6c*EM22 + add122(1, t6c, hu(t6c))
    As = tm(s -> s*(Arb(c1h) + Arb(c1m) + s*s*(Arb(c3) + s*s*(Arb(c5) + s*s*Arb(c7)))) - sin(PI*s), Y, 30.0, 1.0; K = 14)
    Ac = tm(s -> 1 + s*s*(Arb(k2) + s*s*(Arb(k4) + s*s*Arb(k6))) - cos(PI*s), Y, 30.0, 1.0; K = 14)
    (ub(As) + esy, ub(Ac) + ecy, Ym*t6m*(1 + big(2.0)^-40), 1 + t6c)
end
# relative error / allowed eps, for table index idx and quadrant parity q
function piece(idx, q, a::BigFloat, b::BigFloat)
    Y = ival(a, b)
    es, ec, sym, cym = scpiece(Y)
    sah, cah, sam, cam = ST[1, idx+1], ST[2, idx+1], ST[3, idx+1], ST[4, idx+1]
    A = PI*Arb(idx)/128
    esa = ub(Arb(sah) + Arb(sam) - sin(A)); eca = ub(Arb(cah) + Arb(cam) - cos(A))
    SA = abs(sah) + abs(sam); CA = abs(cah) + abs(cam)
    cy = ub(cos(PI*Y)); sy = ub(abs(sin(PI*Y)))
    if iseven(q)       # sa*cy + ca*sy
        tru = sin(A + PI*Y); p1 = SA*cym; p2 = CA*sym
        e = SA*ec + esa*cy + esa*ec + CA*es + eca*sy + eca*es + EM22*(p1 + p2)
    else               # ca*cy - sa*sy
        tru = cos(A + PI*Y); p1 = CA*cym; p2 = SA*sym
        e = CA*ec + eca*cy + eca*ec + SA*es + esa*sy + esa*es + EM22*(p1 + p2)
    end
    m1 = p1*(1 + big(2.0)^-50); m2 = p2*(1 + big(2.0)^-50)
    e += add22cond(m1, hu(m1), m2, hu(m2))
    [e/lb(abs(tru))/ERN]
end

function run()
    @printf("  sinpi_rn: the constant %s allows relative error 2^%.3f\n", RNC, l2(ERN))
    worst = (big(0.0), big(0.0))
    h = big(1)/256
    for q in 0:1, idx in 0:63
        pcs = idx == 0 && q == 0 ? [geopieces(big(2.0)^-960, h, 2); [(-b, -a) for (a, b) in geopieces(big(2.0)^-960, h, 2)]] :
                                   linpieces(-h, h, 32)
        res = bnb((a, b) -> piece(idx, q, a, b), pcs; relw = big(2.0)^-12, maxit = 300)
        worst = supall(worst, res[1])
    end
    [Row("sinpi_rn first step RN relative", worst[1]*ERN, ERN, "$(ctext(F, "rh==rh+1.00001*rm")) constant $RNC, allows")]
end
end
