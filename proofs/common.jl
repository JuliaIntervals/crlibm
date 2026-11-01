using Arblib, Printf

setprecision(Arb, 300)
setprecision(BigFloat, 512)

const SRC = normpath(joinpath(@__DIR__, ".."))

# ---------------------------------------------------------------- C sources

const CACHE = Dict{String,Tuple{String,String}}()

function csource(file)
    get!(CACHE, file) do
        raw = read(joinpath(SRC, file), String)
        s = replace(raw, r"/\*.*?\*/"s => m -> replace(m, r"[^\n]" => " "))
        s = replace(s, r"//[^\n]*" => m -> " "^length(m))
        (raw, s)
    end[2]
end

lineof(file, pos) = 1 + count(==('\n'), csource(file)[1:pos])

function cfind(file, rx)
    m = match(rx, csource(file))
    m === nothing && error("pattern $rx not found in $file")
    m
end

# value of `#define NAME value`, with its location
function cdef(file, name)
    m = cfind(file, Regex("#define\\s+$name\\s+(\\S+)"))
    (m[1], "$file:$(lineof(file, m.offset))")
end
cdouble(file, name) = parse(Float64, cdef(file, name)[1])
cint(file, name) = parse(Int, cdef(file, name)[1])
cloc(file, name) = cdef(file, name)[2]

# value of a local `double NAME = value;` inside a function body
function clocal(file, name)
    m = cfind(file, Regex("double\\s+$name\\s*=\\s*([^;\\s]+)\\s*;"))
    parse(Float64, m[1])
end

const NUM = r"[-+]?(0x[0-9A-Fa-f]+(\.[0-9A-Fa-f]*)?(p[-+]?\d+)?|\d+\.?\d*([eE][-+]?\d+)?|\.\d+([eE][-+]?\d+)?)"

# numeric tokens of the k-th initializer `name[...] = { ... };` in file
function ctokens(file, name; occurrence = 1)
    s = csource(file)
    ms = collect(eachmatch(Regex("\\b$name\\s*(\\[[^\\]]*\\]\\s*)+=\\s*\\{"), s))
    m = ms[occurrence]
    i = m.offset + length(m.match) - 1
    depth = 0; j = i
    for k in i:lastindex(s)
        s[k] == '{' && (depth += 1)
        s[k] == '}' && (depth -= 1)
        if depth == 0; j = k; break; end
    end
    [x.match for x in eachmatch(NUM, s[i:j])]
end
cdoubles(file, name; kw...) = parse.(Float64, ctokens(file, name; kw...))

# db_number tables written as {{0xHI,0xLO}} (big-endian block, first in the file)
dbword(hi, lo) = reinterpret(Float64, (UInt64(hi) << 32) | UInt64(lo))
function cdbnumbers(file, name)
    t = parse.(UInt32, ctokens(file, name; occurrence = 1))
    u = parse.(UInt32, ctokens(file, name; occurrence = 2))
    @assert t[1:2:end] == u[2:2:end] && t[2:2:end] == u[1:2:end] "$name: endian blocks differ"
    [dbword(t[k], t[k+1]) for k in 1:2:length(t)]
end
function cdbnumber(file, name)
    ms = collect(eachmatch(Regex("db_number\\s+const\\s+$name\\s*=\\s*\\{\\{\\s*(0x\\w+)\\s*,\\s*(0x\\w+)\\s*\\}\\}"), csource(file)))
    @assert length(ms) == 2
    h, l = parse(UInt32, ms[1][1]), parse(UInt32, ms[1][2])
    @assert (h, l) == (parse(UInt32, ms[2][2]), parse(UInt32, ms[2][1]))
    dbword(h, l)
end

# line of the first occurrence of a literal piece of text
function ctext(file, text)
    raw = read(joinpath(SRC, file), String)
    i = findfirst(text, raw)
    i === nothing && error("\"$text\" not found in $file")
    "$file:$(1 + count(==('\n'), raw[1:first(i)]))"
end

# ---------------------------------------------------------------- balls and floating-point errors

const A0 = Arb(0)
ub(x::Arb) = BigFloat(Arblib.abs_ubound(Arf, x))
lb(x::Arb) = BigFloat(Arblib.abs_lbound(Arf, x))
rad(r) = Arb((-BigFloat(r), BigFloat(r)))
ival(a, b) = Arb((BigFloat(a), BigFloat(b)))
mid(X::Arb) = Arb(Arblib.midref(X))
radius(X::Arb) = ub(Arb(Arblib.radref(X)))
l2(x) = Float64(log2(BigFloat(x)))

# bound on |RN(v) - v| for every |v| <= M
function hu(M::BigFloat)
    M == 0 && return big(0.0)
    e = exponent(M)
    e < -1022 && return big(2.0)^-1075
    big(2.0)^(e - 53)
end
hu(M::Arb) = hu(ub(M))
hu(M::Real) = hu(BigFloat(M))

# exact value v (ball over the inputs) and bound e on |computed - v|
struct FQ; v::Arb; e::BigFloat; end
FQ(v::Arb) = FQ(v, big(0.0))
FQ(x::Float64) = FQ(Arb(x), big(0.0))
mag(a::FQ) = ub(a.v) + a.e
fround(v, e, exact) = FQ(v, exact ? e : e + hu(ub(v) + e))
fmul(a::FQ, b::FQ; exact = false) = fround(a.v*b.v, ub(a.v)*b.e + ub(b.v)*a.e + a.e*b.e, exact)
fadd(a::FQ, b::FQ; exact = false) = fround(a.v + b.v, a.e + b.e, exact)

# ---------------------------------------------------------------- double-double operators (crlibm_private.h)
# Arguments are bounds on magnitudes; results bound the absolute error of the operator.

function add22(xh, xl, yh, yl)            # Add22, needs |xh| >= |yh|
    m1 = hu(xh + yh) + yl; e1 = hu(m1)
    e1 + hu(m1 + e1 + xl)
end
function add22cond(xh, xl, yh, yl)        # Add22Cond with AVOID_BRANCHES
    e3 = hu(xl + yl)
    e3 + hu(hu(xh + yh) + xl + yl + e3)
end
add122(a, bh, bl) = hu(hu(a + bh) + bl)   # Add122, needs |a| >= |bh|
function mul122(a, bh, bl)
    t3 = a*bl; e3 = hu(t3)
    e3 + hu(hu(a*bh) + t3 + e3)
end
function mul22(xh, xl, yh, yl)
    t1 = xh*yl; t2 = xl*yh; e1 = hu(t1); e2 = hu(t2)
    m3 = t1 + t2 + e1 + e2; e3 = hu(m3)
    xl*yl + e1 + e2 + e3 + hu(hu(xh*yh*(1 + big(2.0)^-52)) + m3 + e3)
end
function muladd22(ch, cl, ah, al, bh, bl) # MulAdd22, needs |ch| >= |ah*bh|
    p = ah*bh*(1 + big(2.0)^-52)
    m7 = hu(p) + cl; e7 = hu(m7)
    m8 = hu(ch + p) + m7 + e7; e8 = hu(m8)
    e5 = hu(ah*bl); e6 = hu(al*bh)
    m9 = ah*bl + al*bh + e5 + e6; e9 = hu(m9)
    al*bl + e5 + e6 + e7 + e8 + e9 + hu(m8 + e8 + m9 + e9)
end
# Div22: x = xh+xl, y = yh+yl, |xh/yh| <= qm, |yh| >= yhlo
function div22(qm, xl, yhlo, yhm, yl)
    r = yhm*hu(qm)
    ch = qm*(1 + big(2.0)^-52)
    m3 = r + xl; e3 = hu(m3)
    mp = ch*yl; ep = hu(mp)
    m4 = m3 + e3 + mp + ep; e4 = hu(m4)
    (m3 + mp)*yl/(yhlo*(yhlo - yl)) + (e3 + ep + e4)/yhlo + hu(m4/yhlo)
end

# relative error of Mul22 for normalized inputs, used below as a constant
const EM22 = big(2.0)^-102.5
function mul22_relative()
    N = 64; w = big(1.0)/N; m = big(0.0)
    for i in 0:N-1, j in 0:N-1
        x0, y0 = 1 + i*w, 1 + j*w
        x1, y1 = min(x0 + w, big(prevfloat(2.0))), min(y0 + w, big(prevfloat(2.0)))
        u = big(2.0)^-53
        m = max(m, mul22(x1, u, y1, u)/(x0*y0*(1 - u)^2))
    end
    m
end
@assert mul22_relative() <= EM22

# ---------------------------------------------------------------- rounding tests

# lower bound of 2^(E(yh)+1) in TEST_AND_RETURN_* when |yh + yl - f| <= e
u53(F, e) = big(2.0)^(exponent((F - e)/(1 + big(2.0)^-52)) + 1)
metrics(F, e) = [e/u53(F, e), e/F]

# compute_rn_constant from maple/common-procedures.mpl
function rn_constant(eps::Arb)
    k = floor(Int, -log2(BigFloat(Arblib.midref(eps)))) - 53
    (1 + Arb(2)^54*eps/(1 - eps - Arb(2)^(1 - k)))/(1 - Arb(2)^-53)
end
# largest eps whose constant does not exceed e (checked in ball arithmetic)
function eps_rn(e)
    A = big(e)*(1 - big(2.0)^-53) - 1
    ep = A/(big(2.0)^54 + A)
    for _ in 1:5
        k = floor(Int, -log2(ep)) - 53
        ep = A*(1 - big(2.0)^(1 - k))/(big(2.0)^54 + A)
    end
    ep *= 1 - big(2.0)^-200
    c = rn_constant(Arb(ep))
    @assert BigFloat(Arblib.ubound(Arf, c)) <= big(e) "eps_rn($e) not verified"
    ep
end

# ---------------------------------------------------------------- approximation errors

# enclosure of f over X by a Taylor model of degree K at the centre of X;
# M bounds |f| on the complex disc of radius rho around the centre
function tm(f, X::Arb, M, rho; K = 16)
    c = mid(X); r = radius(X); @assert r < rho/4
    D = f(ArbSeries((c, 1), degree = K))
    t = rad(r); v = A0
    for k in K:-1:0; v = v*t + D[k]; end
    v + rad(big(M)*(r/big(rho))^(K + 1)/(1 - r/big(rho)))
end

# ---------------------------------------------------------------- covers and branch and bound

geopieces(lo, hi, ratio) = (p = Tuple{BigFloat,BigFloat}[]; a = big(lo);
    while a < hi; b = min(big(hi), a*ratio); push!(p, (a, b)); a = b; end; p)
linpieces(lo, hi, n) = [(big(lo) + (big(hi) - lo)*k/n, big(lo) + (big(hi) - lo)*(k + 1)/n) for k in 0:n-1]

# the pieces, sorted, cover [lo, hi] without gaps
function covers(pieces, lo, hi)
    p = sort(pieces)
    p[1][1] <= lo && p[end][2] >= hi && all(p[k][2] >= p[k+1][1] for k in 1:length(p)-1)
end

# f(a, b) -> vector of upper bounds over [a, b]; returns, per objective,
# (sup bound over the cover, centre of the worst piece). Refining keeps a cover.
function bnb(f, pieces; relw = big(2.0)^-18, maxit = 4000)
    cache = Dict{Tuple{BigFloat,BigFloat},Vector{BigFloat}}()
    F(a, b) = get!(() -> f(a, b), cache, (a, b))
    nobj = length(F(pieces[1]...))
    res = Tuple{BigFloat,BigFloat}[]
    for o in 1:nobj
        leaves = [(a, b, F(a, b)[o]) for (a, b) in pieces]
        it = 0
        while true
            k = argmax([l[3] for l in leaves]); a, b, w = leaves[k]
            if (b - a) <= relw*max(abs(a), abs(b)) || it >= maxit
                push!(res, (w, (a + b)/2)); break
            end
            m = (a + b)/2; deleteat!(leaves, k)
            push!(leaves, (a, m, F(a, m)[o]), (m, b, F(m, b)[o])); it += 1
        end
    end
    res
end
supall(rs...) = reduce((x, y) -> x[1] >= y[1] ? x : y, rs)

# ---------------------------------------------------------------- results

struct Row
    name::String      # what is bounded
    bound::BigFloat   # proven bound
    cst::BigFloat     # bound the code relies on
    src::String       # where the code takes it from
end
rowdir(name, b, file, cname) = Row(name, b, big(cdouble(file, cname)), "$(cloc(file, cname)) $cname")
rowrn(name, b, file, cname) = (c = cdouble(file, cname);
    Row(name, b, eps_rn(c), "$(cloc(file, cname)) $cname = $(c), allows"))
