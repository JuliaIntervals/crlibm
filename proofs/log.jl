# log.c quick phase: log, log2, log10 (all rounding modes)

module LogProof
include("common.jl")

const H = "log.h"
const T = reshape(cdoubles(H, "argredtable"), 4, 256)
const RI = [Float64(Float32(T[1, i])) for i in 1:256]   # ri is a float
const LOGIH, LOGIM = T[2, :], T[3, :]
const C = [cdouble(H, "p_coeff_$(k)h") for k in 1:6]
const log2h, log2m = cdouble(H, "log2h"), cdouble(H, "log2m")
const log210h, log210m = cdouble(H, "log210h"), cdouble(H, "log210m")
const RH2, RL2 = cdouble(H, "RECPRLOG2H"), cdouble(H, "RECPRLOG2L")
const RH10, RM10 = cdouble(H, "RECPRLOG10H"), cdouble(H, "RECPRLOG10M")
@assert C[1] == 1 && C[2] == -0.5
@assert cint(H, "L") == 8 && cint(H, "MAXINDEX") == 106 && cint(H, "INDEXMASK") == 255

sigbits(x) = (x == 0 ? 0 : 53 - trailing_zeros(reinterpret(UInt64, abs(x)) | (UInt64(1) << 52)))
# yh*ri exact (yh has 21 bits) and yl*ri exact (yl has at most 32 bits)
@assert all(sigbits(r) <= 21 for r in RI)
# E*log2h, E*log2m, E*log210h, E*log210m exact for |E| <= 1075 (11 bits)
@assert all(sigbits(c) <= 42 for c in (log2h, log2m, log210h, log210m))
@assert all(abs(l) <= hu(abs(h)) for (h, l) in ((RH2, RL2), (RH10, RM10)))

const LN2, LN10 = log(Arb(2)), log(Arb(10))

# e(z) = z + c2 z^2 + ... + c6 z^6 - log1p(z), as a series
const QK = Dict(k => (k <= 6 ? Arb(C[k]) - Arb((-1)^(k + 1))/k : Arb((-1)^k)/k) for k in 3:60)
function ea(z::Arb)
    s = A0
    for k in 60:-1:3; s = s*z + QK[k]; end
    zm = ub(z)
    s*z*z*z + rad(zm^61/(61*(1 - zm)))
end
function dea(z::Arb)
    s = A0
    for k in 60:-1:3; s = s*z + k*QK[k]; end
    zm = ub(z)
    s*z*z + rad(zm^60/(1 - zm)^2)
end
eaC(Z::Arb) = ea(mid(Z)) + dea(Z)*rad(radius(Z))

# error of logTabPolyh + logTabPolyl against log(y) for zh in Z (table index i)
# returns (signed part, bound of the rest, log(y), bound of |logTabPoly|)
function ypart(i, Z::Arb)
    j = (i & 255) + 1; ri, logih, logim = RI[j], LOGIH[j], LOGIM[j]
    Zm = ub(Z)
    Hl = ri == 1.0 ? big(0.0) : hu(Zm)            # |zl|; zl = 0 when ri = 1
    zh = FQ(Z)
    s = fmul(zh, zh)
    p35 = fadd(FQ(C[3]), fmul(s, FQ(C[5])))
    p46 = fadd(FQ(C[4]), fmul(s, FQ(C[6])))
    zc = fmul(s, zh)
    zsh = fmul(FQ(C[2]), s; exact = true)
    z4 = fmul(s, s)
    p36 = fadd(fmul(zc, p35), fmul(z4, p46))
    t = FQ(zsh.v, zsh.e + hu(ub(zsh.v) + zsh.e + Hl))   # + zl, zl kept apart
    pU = FQ(t.v + p36.v, t.e + p36.e + hu(ub(t.v + p36.v) + Hl + t.e + p36.e))
    # Add12(ph, pl, zh, pUpper): |pUpper| <= |zh|*(|zh|/2 + 2^-53 + zh^2)*(1+2^-53)^3 < |zh|
    @assert (Zm*(1 + big(2.0)^-53)/2 + big(2.0)^-53 + Zm^2)*(1 + big(2.0)^-53)^3 < 1
    det = Arb(logih) + Arb(logim) + log(Arb(ri)) + eaC(Z)
    zlt = Hl*Zm/(1 - Zm) + 2*Hl^2
    phm = Zm + ub(pU.v) + Hl + pU.e
    if logih == 0.0
        d1 = big(0.0)
    else
        @assert abs(logih) >= phm*(1 + big(2.0)^-50)
        d1 = add22(abs(big(logih)), abs(big(logim)), phm, hu(phm))
    end
    f = log1p(Z + rad(Hl)) - log(Arb(ri))
    R = Arb(logih) + Z + pU.v
    return det, zlt + pU.e + d1, f, ub(R) + Hl + pU.e + big(2.0)^-100
end

# E = 0: (directed metric, relative error) for log, log2, log10
function ratios(i, Z)
    det, rest, f, Rm = ypart(i, Z)
    errR = ub(det) + rest
    out = BigFloat[]
    for fn in 0:2
        if fn == 0
            err = errR; F = f
        else
            XH, XL, L = fn == 1 ? (RH2, RL2, LN2) : (RH10, RM10, LN10)
            rho = ub(Arb(XH) + Arb(XL) - 1/L)
            err = (abs(XH) + abs(XL))*errR + ub(f)*rho + mul22(abs(big(XH)), abs(big(XL)), Rm, hu(Rm))
            F = f/L
        end
        Fl = lb(F)
        Fl <= err && (append!(out, [big(Inf), big(Inf)]); continue)
        append!(out, metrics(Fl, err))
    end
    out
end

# zh range of table cell i (index before masking, 0..256)
function cellrange(i)
    lo = i*4096 - 2048; hi = i*4096 + 2048
    i == 0 && (lo = 0); i == 256 && (hi = 2^20)
    ylo = 1 + big(lo)/2^20; yhi = 1 + big(hi)/2^20 - big(2.0)^-52
    if i >= 106; ylo /= 2; yhi /= 2; end
    ri = big(RI[(i & 255) + 1])
    (ylo*ri - 1, yhi*ri - 1)
end

function pieces_of(i)
    zlo, zhi = cellrange(i)
    p = Tuple{BigFloat,BigFloat}[]
    for (a, b) in ((zlo, min(zhi, big(0.0))), (max(zlo, big(0.0)), zhi))
        a >= b && continue
        if a == 0 || b == 0
            s = a == 0 ? 1 : -1; m = max(abs(a), abs(b)); k = exponent(m)
            push!(p, s > 0 ? (big(2.0)^k, m) : (-m, -big(2.0)^k))
            for j in k-1:-1:-56
                push!(p, s > 0 ? (big(2.0)^j, big(2.0)^(j+1)) : (-big(2.0)^(j+1), -big(2.0)^j))
            end
        else
            append!(p, linpieces(a, b, 16))
        end
    end
    # with ri = 1, zh = 0 only for x = 2^E; E = 0 gives x = 1, where every step is exact
    RI[(i & 255) + 1] != 1.0 && zlo < 0 < zhi && push!(p, (-big(2.0)^-56, big(2.0)^-56))
    @assert covers(p, i == 0 ? big(2.0)^-56 : zlo, zhi)
    p
end

# E != 0: the y-part error is bounded once per cell, then combined with E*log(2)
function cellsup(i)
    zlo, zhi = cellrange(i)
    eY = big(0.0); fM = big(0.0); RM = big(0.0)
    for (a, b) in linpieces(zlo, zhi, 64)
        det, rest, f, Rm = ypart(i, ival(a, b))
        eY = max(eY, ub(det) + rest); fM = max(fM, ub(f)); RM = max(RM, Rm)
    end
    eY, fM, RM
end

function Enonzero()
    rep2 = ub(Arb(log2h) + Arb(log2m) - LN2)
    rep10 = ub(Arb(log210h) + Arb(log210m) - LN2/LN10)
    rho2 = ub(Arb(RH2) + Arb(RL2) - 1/LN2); rho10 = ub(Arb(RH10) + Arb(RM10) - 1/LN10)
    w = zeros(BigFloat, 6)
    for i in 0:256
        eY, fM, RM = cellsup(i)
        for E in -1074:1024
            E == 0 && continue
            aE = abs(big(E))
            xh = aE*big(log2h); @assert xh >= RM
            err = eY + aE*rep2 + add22(xh, hu(xh), RM, hu(RM))
            r1 = metrics(lb(aE*LN2) - fM, err)
            R2 = (abs(RH2) + abs(RL2))*RM; @assert aE >= R2
            err2 = (abs(RH2) + abs(RL2))*eY + fM*rho2 + mul22(abs(big(RH2)), abs(big(RL2)), RM, hu(RM))
            err2 += add122(aE, R2 + err2, hu(R2 + err2))
            r2 = metrics(aE - ub(fM/LN2), err2)
            R10 = (abs(RH10) + abs(RM10))*RM
            xh = aE*big(log210h); @assert xh >= R10
            err10 = (abs(RH10) + abs(RM10))*eY + fM*rho10 + mul22(abs(big(RH10)), abs(big(RM10)), RM, hu(RM))
            err10 += aE*rep10 + add22(xh, hu(xh), R10 + err10, hu(R10 + err10))
            r10 = metrics(lb(aE*LN2/LN10) - ub(fM/LN10), err10)
            w = max.(w, [r1; r2; r10])
        end
    end
    w
end

function run()
    best = [(big(0.0), big(0.0)) for _ in 1:6]
    for i in 0:256
        res = bnb((a, b) -> ratios(i, ival(a, b)), pieces_of(i); relw = big(2.0)^-16, maxit = 3000)
        best = [supall(best[o], res[o]) for o in 1:6]
    end
    wE = Enonzero()
    b = [max(best[o][1], wE[o]) for o in 1:6]
    @printf("  log.c: worst bounds for E != 0: 2^%.2f (relative 2^%.2f)\n", l2(maximum(wE[1:2:5])), l2(maximum(wE[2:2:6])))
    [rowdir("log  directed", b[1], H, "RDROUNDCST"),
     rowdir("log2 directed", b[3], H, "RDROUNDCST"),
     rowdir("log10 directed", b[5], H, "RDROUNDCST"),
     rowrn("log  RN relative", b[2], H, "RNROUNDCST"),
     rowrn("log2 RN relative", b[4], H, "RNROUNDCST"),
     rowrn("log10 RN relative", b[6], H, "RNROUNDCST")]
end
end
