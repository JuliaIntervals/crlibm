# Runs every proof and checks each proven bound against the bound the C code relies on.
# Usage: julia --project=proofs proofs/run.jl

using Printf

const FILES = ["log.jl" => :LogProof, "asincos.jl" => :AsinProof, "trigo.jl" => :TrigProof, "atan.jl" => :AtanProof,
               "expm1.jl" => :Expm1Proof, "exp.jl" => :ExpProof, "sinpi.jl" => :SinpiProof]

f2(x) = Float64(log2(BigFloat(x)))
rev = try readchomp(`git -C $(joinpath(@__DIR__, "..")) rev-parse --short HEAD`) catch; "unknown" end
println("crlibm first-step error bounds, sources at $(normpath(joinpath(@__DIR__, ".."))) (commit $rev)")
failed = 0
t0 = time()
for (file, mod) in FILES
    t = time()
    println("\n", file)
    rows = try
        include(joinpath(@__DIR__, file))
        Base.invokelatest(getfield(getfield(Main, mod), :run))
    catch err
        println("  FAILED: ", sprint(showerror, err))
        global failed += 1
        continue
    end
    for r in rows
        ok = r.bound < r.cst
        ok || (global failed += 1)
        @printf("  %-42s 2^%8.3f %s 2^%8.3f  margin %6.3f bit  %s  %s\n",
                r.name, f2(r.bound), ok ? " <" : ">=", f2(r.cst), f2(r.cst) - f2(r.bound), ok ? "ok  " : "FAIL", r.src)
    end
    @printf("  (%.0f s)\n", time() - t)
end
@printf("\ntotal %.0f s; %s\n", time() - t0, failed == 0 ? "all bounds proven" : "$failed FAILURES")
exit(failed == 0 ? 0 : 1)
