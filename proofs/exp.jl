# exp-td.c quick phase, with the Cody-Waite reduction Add12Cond(rh, rm, x - kd*Log2h, -kd*Log2l)

module ExpProof
include("common.jl")

const H = "exp-td.h"
const T1 = reshape(cdoubles(H, "twoPowerIndex1"), 3, 64); const T2 = reshape(cdoubles(H, "twoPowerIndex2"), 3, 64)
const CINV = cdouble(H, "log2InvMult2L")
const C3, C4 = cdouble(H, "c3"), cdouble(H, "c4")
const OVF, DEN = cdouble(H, "OVRFLWBOUND"), cdouble(H, "DENORMBOUND")
const Log2h, Log2l = clocal("exp-td.c", "Log2h"), clocal("exp-td.c", "Log2l")
@assert cint(H, "L") == 12 && cint(H, "INDEXMASK1") == 0x3f && cint(H, "INDEXMASK2") == 0xfc0
const LN2 = log(Arb(2))

sigbits(x) = 53 - trailing_zeros(reinterpret(UInt64, abs(x)) | (UInt64(1) << 52))
const XMAX = max(abs(big(OVF)), abs(big(DEN)))      # quick phase only for DENORMBOUND < x <= OVRFLWBOUND
const KMAX = XMAX*CINV + 1
const REP = ub(Arb(Log2h) + Arb(Log2l) - LN2/4096)
# kd*Log2h exact, RN(kd*Log2l) rounded once, Log2h + Log2l against ln(2)/4096
const ERED = hu(KMAX*abs(Log2l)) + KMAX*REP
const RMAX = ((0.5 + hu(KMAX) + XMAX*ub(Arb(CINV) - 4096/LN2))*ub(LN2/4096) + ERED)*(1 + big(2.0)^-40)

const ETAB = let
    e1 = maximum(ub((Arb(T1[1, i+1]) + Arb(T1[2, i+1]))/Arb(2)^(Arb(i)/4096) - 1) for i in 0:63)
    e2 = maximum(ub((Arb(T2[1, i+1]) + Arb(T2[2, i+1]))/Arb(2)^(Arb(i)/64) - 1) for i in 0:63)
    (1 + e1)*(1 + e2)*(1 + EM22) - 1
end

function quickR(R::Arb)
    Rm = ub(R); Hm = hu(Rm)
    rh = FQ(R)
    sq = fmul(rh, rh); c3 = fmul(FQ(C3), rh); sqh = fmul(FQ(0.5), sq; exact = true)
    hp = fadd(fmul(c3, sq), fmul(FQ(C4), fmul(sq, sq))); hpws = fadd(sqh, hp)
    t8e = hpws.e + hu(ub(hpws.v) + hpws.e + Hm)
    t9m = ub(R + hpws.v) + Hm + t8e
    t9e = t8e + hu(t9m)
    A = tm(s -> exp(s) - 1 - (s + s*s/2 + Arb(C3)*s*s*s + Arb(C4)*s*s*s*s), R, 4.0, 1.0; K = 14)
    ETAB*(1 + t9m) + ub(A) + Hm*Rm*1.001 + Hm^2 + t9e + 2*big(2.0)^-53*t9m + big(2.0)^-104 + ERED*1.0002
end

function run()
    @assert sigbits(Log2h) + ceil(Int, log2(KMAX)) <= 53
    # x - kd*Log2h exact: for kd != 0 both are multiples of 2^(emin-52), emin >= exponent(Log2h/2 - dmax),
    # and |x - kd*Log2h| <= dmax < 2^(emin+1); for kd = 0 it is x
    dmax = (0.5 + hu(KMAX) + XMAX*ub(Arb(CINV) - 4096/LN2))*ub(LN2/4096) + KMAX*ub(Arb(Log2h) - LN2/4096)
    @assert dmax < big(2.0)^(exponent(Log2h - dmax) + 1)
    c = (big(1.0) + RMAX)
    @assert abs(1/big(6) - C3)*c^3 + abs(1/big(24) - C4)*c^4 + ub(exp(Arb(c))) - sum(c^k/factorial(big(k)) for k in 0:4) <= 4
    m = maximum(quickR(ival(a, b))/lb(exp(ival(a, b))) for (a, b) in linpieces(-RMAX, RMAX, 512))
    @printf("  exp: table 2^%.2f, |r| <= %.6e, exp quick relative error <= 2^%.3f\n", l2(ETAB), Float64(RMAX), l2(m))
    # polyTblh + polyTblm is in [1/2, 4): |error|/2^(E(polyTblh)+1) <= relative error * (1 + 2^-50)
    [Row("exp reduction |r - (rh + rm)|", ERED, big(2.0)^-77, "$(ctext("exp-td.c", "accurate to 2^-77")) comment"),
     rowdir("exp directed", m*(1 + big(2.0)^-50), H, "RDROUNDCST"),
     rowrn("exp RN relative", m, H, "ROUNDCST")]
end
end
