# ===========================================================================
#  yl_collapse_stages23.jl — SELF-CONTAINED. Stages 2 and 3 of the BALFE-M -> LFE-M
#  operator collapse against Yang & Liu (2024), JFM 999 A32 + supplementary material.
#
#      julia --project=. test/yl_collapse_wip/yl_collapse_stages23.jl
#
#  Stage 1 (vertical velocity vs supplementary A) is the committed gate,
#  test/test_yl_collapse.jl.
#
#  STAGE 2  our (L,θ)+(N,Θ) package  ==  Dw/Dt        -> verifies all 8 Θ_kj
#  STAGE 3  our TENSOR-FORM residual ==  H ∫φᵢ R dσ   -> verifies the SOLVER's tensors
#
#  ⚠ Both run on a CONTINUITY-CONSISTENT state. pDerivation.tex eliminates ∂H/∂t via
#  the depth-integrated continuity equation, so our package is an identity ONLY on
#  that constraint; Yang & Liu keep H^(0,1) explicit and are state-independent. A
#  state with H and u_j prescribed independently makes OUR side fail by O(eps^2) for
#  reasons that are an artefact of the test, not a defect. See README.md.
# ===========================================================================
using ForwardDiff, Printf, LinearAlgebra
using GridapBALFEM                      # assemble_vertical_tensors (stage 3)

const d0 = 3.5
const AMP  = Ref(1.0)                   # wave-amplitude scale (bathymetry NOT scaled)
const FLATBED = Ref(false)

hbed(x)    = FLATBED[] ? d0 : d0 - 0.8*exp(-((x-12.0)/4.0)^2)

# ---------------- p=1 basis, explicit (not taken from the solver) ----------
function phi(c, m, s)
    n = length(c)
    (m > 1 && c[m-1] <= s <= c[m]) && return (s - c[m-1])/(c[m] - c[m-1])
    (m < n && c[m] <= s <= c[m+1]) && return (c[m+1] - s)/(c[m+1] - c[m])
    return zero(s)
end
function dphi(c, m, s)              # φ'_m : piecewise constant
    n = length(c)
    (m > 1 && c[m-1] <= s <= c[m]) && return 1/(c[m] - c[m-1])
    (m < n && c[m] <= s <= c[m+1]) && return -1/(c[m+1] - c[m])
    return 0.0
end
const G8, W8 = let n=8
    i=1:(n-1); b=@. i/sqrt(4.0*i^2-1.0); J=diagm(-1=>collect(b),1=>collect(b))
    F=eigen(Symmetric(J)); (F.values, 2.0.*(F.vectors[1,:].^2))
end
"∫_a^b f, composite 8-pt Gauss on the mesh sub-intervals (exact for our degrees)"
function gquad(f, c, a, b)
    a >= b && return 0.0
    pts = sort(unique(vcat(a, b, [x for x in c if a < x < b])))
    tot = 0.0
    for i in 1:length(pts)-1
        lo, hi = pts[i], pts[i+1]; hw = (hi-lo)/2; mid = (lo+hi)/2
        for (g,w) in zip(G8,W8); tot += hw*w*f(mid+hw*g); end
    end
    tot
end
phi_int(c,m,s) = gquad(t->phi(c,m,t), c, 0.0, s)
Phi_w(c,m)     = gquad(t->phi(c,m,t), c, 0.0, 1.0)

# ---------------- prescribed state, now a function of (x,t) ----------------
const d0 = 3.5
const AMP = Ref(1.0)
const FLATBED = Ref(false)
hbed(x)     = FLATBED[] ? d0 : d0 - 0.8*exp(-((x-12.0)/4.0)^2)
etaf(x,t)   = AMP[]*(0.09*cos(1.1x + 0.3 - 0.8t) + 0.03*sin(2.3x + 0.5t))
Hf(x,t)     = hbed(x) + etaf(x,t)
uf(x,t,m)   = AMP[]*(0.4 + 0.15m)*sin((0.7 + 0.13m)*x + 0.4m - (0.6 + 0.05m)*t)

dx(f,x,t)   = ForwardDiff.derivative(z->f(z,t), x)
dt_(f,x,t)  = ForwardDiff.derivative(z->f(x,z), t)
dxx(f,x,t)  = ForwardDiff.derivative(z->ForwardDiff.derivative(y->f(y,t), z), x)

# ---------------- THEIR w coefficients, as functions of (x,t) --------------
function wcoef(c, x, t)                      # returns w0,w1,w2 (vectors over elements)
    M = length(c)-1
    H  = Hf(x,t); Hx = dx(Hf,x,t); hx = ForwardDiff.derivative(hbed,x)
    u  = [uf(x,t,m) for m in 1:M+1]
    ux = [dx((z,τ)->uf(z,τ,m), x, t) for m in 1:M+1]
    w1 = [ (hx*u[k+1] - hx*u[k] + c[k+1]*H*ux[k] - c[k]*H*ux[k+1])/(c[k]-c[k+1]) for k in 1:M ]
    w2 = [ (Hx*u[k] - Hx*u[k+1] + H*ux[k+1] - H*ux[k])/(2*(c[k]-c[k+1]))          for k in 1:M ]
    w0 = zeros(typeof(H), M); w0[1] = -u[1]*hx
    for k in 2:M
        s = -u[1]*hx
        for m in 1:k-1
            s += c[m+1]*(w1[m]-w1[m+1]) + c[m+1]^2*(w2[m]-w2[m+1])
        end
        w0[k] = s
    end
    return w0, w1, w2
end
wn(c,x,t,k,n) = (W=wcoef(c,x,t); n==0 ? W[1][k] : n==1 ? W[2][k] : W[3][k])


const T0   = 0.0
const CMESH = Ref([0.0,0.728,1.0])

H0f(x)      = hbed(x) + AMP[]*(0.09*cos(1.1x + 0.3) + 0.03*sin(2.3x))
u0f(x,t,m)  = AMP[]*(0.4 + 0.15m)*sin((0.7 + 0.13m)*x + 0.4m - (0.6 + 0.05m)*t)
function Gf(x)
    c = CMESH[]
    -sum(ForwardDiff.derivative(z -> H0f(z)*u0f(z,T0,m), x)*Phi_w(c,m) for m in 1:length(c))
end
Hc(x,t) = H0f(x) + (t - T0)*Gf(x)

#  w, ω and the two material-derivative forms, all on the consistent state
function w_c(c,x,t,s)
    H = Hc(x,t); Hx = ForwardDiff.derivative(z->Hc(z,t),x)
    hx = ForwardDiff.derivative(hbed,x)
    acc = zero(H)
    for m in 1:length(c)
        u = u0f(x,t,m); ux = ForwardDiff.derivative(z->u0f(z,t,m),x)
        acc += u*hx*phi(c,m,s) - u*Hx*s*phi(c,m,s) + (Hx*u + H*ux)*phi_int(c,m,s)
    end
    -acc
end
omega_c(c,x,t,s) = sum(
    (ForwardDiff.derivative(z->Hc(z,t)*u0f(z,t,m), x))*(s*Phi_w(c,m) - phi_int(c,m,s))
    for m in 1:length(c))

"YL form: w_t + w_σ σ_t + u(w_x + w_σ σ_x) + w w_σ σ_z"
function DwDt_YL(c,x,t,s)
    H  = Hc(x,t); Hx = ForwardDiff.derivative(z->Hc(z,t),x)
    Ht = ForwardDiff.derivative(z->Hc(x,z),t); hx = ForwardDiff.derivative(hbed,x)
    u  = sum(u0f(x,t,m)*phi(c,m,s) for m in 1:length(c))
    w  = w_c(c,x,t,s)
    wx = ForwardDiff.derivative(z->w_c(c,z,t,s), x)
    wt = ForwardDiff.derivative(z->w_c(c,x,z,s), t)
    ws = ForwardDiff.derivative(z->w_c(c,x,t,z), s)
    wt + ws*(-(s/H)*Ht) + u*(wx + ws*(hx - s*Hx)/H) + w*ws/H
end
"our form: w_t + u·∇w + (ω/H) w_σ"
function DwDt_ours(c,x,t,s)
    H  = Hc(x,t)
    u  = sum(u0f(x,t,m)*phi(c,m,s) for m in 1:length(c))
    wx = ForwardDiff.derivative(z->w_c(c,z,t,s), x)
    wt = ForwardDiff.derivative(z->w_c(c,x,z,s), t)
    ws = ForwardDiff.derivative(z->w_c(c,x,t,z), s)
    wt + u*wx + (omega_c(c,x,t,s)/H)*ws
end


D1(f,x,t) = ForwardDiff.derivative(z->f(z,t), x)
D2(f,x,t) = ForwardDiff.derivative(z->D1(f,z,t), x)

#  ⚠ The earlier version used a helper d1(f) that returned a NUMBER, so the lambda
#  (z,tau) -> d1(Hu(k)) * u_j(z,tau) froze dx(Hu_k) as a constant and dropped the
#  u_j * dxx(Hu_k) half of N_2 -- a value captured where a function was meant.
function package_fixed(c,x,t,s)
    M1 = length(c)
    H  = Hc(x,t); Hx = D1(Hc,x,t); hx = ForwardDiff.derivative(hbed,x)
    u  = [u0f(x,t,m) for m in 1:M1]
    ut = [ForwardDiff.derivative(z->u0f(x,z,m), t) for m in 1:M1]
    Hu(m)  = (z,tau)->Hc(z,tau)*u0f(z,tau,m)
    Hux=[D1(Hu(m),x,t) for m in 1:M1]; Huxx=[D2(Hu(m),x,t) for m in 1:M1]
    Hut=[D1((z,tau)->Hc(z,tau)*ForwardDiff.derivative(b->u0f(z,b,m),tau),x,t) for m in 1:M1]
    ujhx=[D1((z,tau)->u0f(z,tau,m)*ForwardDiff.derivative(hbed,z),x,t) for m in 1:M1]
    ujHx=[D1((z,tau)->u0f(z,tau,m)*D1(Hc,z,tau),x,t) for m in 1:M1]
    Ph=[Phi_w(c,m) for m in 1:M1]
    tot = 0.0
    for m in 1:M1
        tot += (-ut[m]*hx)*phi(c,m,s) + (ut[m]*Hx)*(s*phi(c,m,s)) + (-Hut[m])*phi_int(c,m,s)
    end
    for k in 1:M1, j in 1:M1
        N=(-u[j]*Huxx[k], D1((z,tau)->D1(Hu(k),z,tau)*u0f(z,tau,j),x,t), -u[k]*ujhx[j],
            u[k]*ujHx[j], -u[k]*Huxx[j], -Hux[j]*u[k]*hx/H, Hux[j]*u[k]*Hx/H, -Hux[j]*Hux[k]/H)
        Th=(s*Ph[k]*phi(c,j,s), Ph[k]*phi_int(c,j,s), phi(c,j,s)*phi(c,k,s),
            s*phi(c,j,s)*phi(c,k,s), phi_int(c,j,s)*phi(c,k,s),
            s*Ph[j]*dphi(c,k,s)-phi_int(c,j,s)*dphi(c,k,s),
            s*Ph[j]*phi(c,k,s)+s^2*Ph[j]*dphi(c,k,s)-phi_int(c,j,s)*phi(c,k,s)-s*phi_int(c,j,s)*dphi(c,k,s),
            s*Ph[j]*phi(c,k,s)-phi_int(c,j,s)*phi(c,k,s))
        for q in 1:8; tot += N[q]*Th[q]; end
    end
    tot
end

function run_stage2()
    println("STAGE 2 - (L,theta)+(N,Theta) package vs Dw/Dt, continuity-consistent state
")
    @printf("%-12s %-6s %-6s %14s %12s
","mesh","amp","x","max|pkg-DwDt|","relative")
    println("-"^58)
    worst = 0.0
    for cc in ([0.0,0.728,1.0],[0.0,0.5,1.0],[0.0,0.4,0.8,1.0],[0.0,0.3,0.6,0.85,1.0])
        CMESH[] = cc
        for e in (1.0, 0.1), x in (3.7, 12.0, 19.4)
            AMP[] = e
            r = 0.0; sc = 0.0
            for s in range(0.03,0.97;length=31)
                k=searchsortedlast(cc,s); (s-cc[k]<1e-3||cc[k+1]-s<1e-3) && continue
                tgt = DwDt_ours(cc,x,T0,s)
                r = max(r, abs(package_fixed(cc,x,T0,s)-tgt)); sc = max(sc, abs(tgt))
            end
            worst = max(worst, r/max(sc,1e-12))
            @printf("M=%d nodes    %-6.2f %-6.1f %14.3e %12.3e
", length(cc),e,x,r,r/max(sc,1e-12))
        end
    end
    AMP[]=1.0
    @printf("
STAGE 2 WORST RELATIVE DIFFERENCE: %.3e

", worst)
end

uh(c,x,t,s) = sum(u0f(x,t,m)*phi(c,m,s) for m in 1:length(c))
etaf_c(x,t) = Hc(x,t) - hbed(x)

"our p_nh/ρ on the consistent state (L·Π + N·Λ), σ-resolved"
function pnh_c(c,x,t,s)
    M1 = length(c)
    H  = Hc(x,t); Hx = D1(Hc,x,t); hx = ForwardDiff.derivative(hbed,x)
    u  = [u0f(x,t,m) for m in 1:M1]
    ut = [ForwardDiff.derivative(z->u0f(x,z,m), t) for m in 1:M1]
    Hu(m)  = (z,τ)->Hc(z,τ)*u0f(z,τ,m)
    Hux=[D1(Hu(m),x,t) for m in 1:M1]; Huxx=[D2(Hu(m),x,t) for m in 1:M1]
    Hut=[D1((z,τ)->Hc(z,τ)*ForwardDiff.derivative(b->u0f(z,b,m),τ),x,t) for m in 1:M1]
    ujhx=[D1((z,τ)->u0f(z,τ,m)*ForwardDiff.derivative(hbed,z),x,t) for m in 1:M1]
    ujHx=[D1((z,τ)->u0f(z,τ,m)*D1(Hc,z,τ),x,t) for m in 1:M1]
    Ph=[Phi_w(c,m) for m in 1:M1]
    acc = zero(H)
    for m in 1:M1
        L=(-ut[m]*hx, ut[m]*Hx, -Hut[m])
        Pi=(gquad(z->phi(c,m,z),c,s,1.0), gquad(z->z*phi(c,m,z),c,s,1.0),
            gquad(z->phi_int(c,m,z),c,s,1.0))
        acc += sum(L[q]*Pi[q] for q in 1:3)
    end
    for k in 1:M1, j in 1:M1
        N2 = D1((z,τ)->D1(Hu(k),z,τ)*u0f(z,τ,j), x, t)
        N=(-u[j]*Huxx[k], N2, -u[k]*ujhx[j], u[k]*ujHx[j], -u[k]*Huxx[j],
           -Hux[j]*u[k]*hx/H, Hux[j]*u[k]*Hx/H, -Hux[j]*Hux[k]/H)
        Th=(z->z*Ph[k]*phi(c,j,z), z->Ph[k]*phi_int(c,j,z), z->phi(c,j,z)*phi(c,k,z),
            z->z*phi(c,j,z)*phi(c,k,z), z->phi_int(c,j,z)*phi(c,k,z),
            z->z*Ph[j]*dphi(c,k,z)-phi_int(c,j,z)*dphi(c,k,z),
            z->z*Ph[j]*phi(c,k,z)+z^2*Ph[j]*dphi(c,k,z)-phi_int(c,j,z)*phi(c,k,z)-z*phi_int(c,j,z)*dphi(c,k,z),
            z->z*Ph[j]*phi(c,k,z)-phi_int(c,j,z)*phi(c,k,z))
        for q in 1:8; acc += N[q]*gquad(Th[q],c,s,1.0); end
    end
    H*acc
end

"pointwise momentum residual R(σ) of (2.12), per unit density"
function Rpoint(c,x,t,s)
    H=Hc(x,t); Hx=D1(Hc,x,t); Ht=ForwardDiff.derivative(z->Hc(x,z),t)
    hx=ForwardDiff.derivative(hbed,x)
    sig_t=-(s/H)*Ht; sig_x=(hx-s*Hx)/H; sig_z=1/H
    u  = uh(c,x,t,s)
    ut = ForwardDiff.derivative(z->uh(c,x,z,s), t)
    ux = ForwardDiff.derivative(z->uh(c,z,t,s), x)
    us = ForwardDiff.derivative(z->uh(c,x,t,z), s)
    w  = w_c(c,x,t,s)
    etax = ForwardDiff.derivative(z->etaf_c(z,t), x)
    px  = ForwardDiff.derivative(z->pnh_c(c,z,t,s), x)
    ps  = ForwardDiff.derivative(z->pnh_c(c,x,t,z), s)
    ut + us*sig_t + u*(ux + us*sig_x) + w*us*sig_z + 9.81*etax + px + ps*sig_x
end

"our tensor-form residual at DOF i, from the SOLVER's tensors"
function Res_ours(V, c, x, t, i)
    M1 = length(c); H = Hc(x,t); Hx = D1(Hc,x,t); hx = ForwardDiff.derivative(hbed,x)
    u  = [u0f(x,t,m) for m in 1:M1]
    ux = [D1((z,τ)->u0f(z,τ,m),x,t) for m in 1:M1]
    ut = [ForwardDiff.derivative(z->u0f(x,z,m), t) for m in 1:M1]
    Hu(m)=(z,τ)->Hc(z,τ)*u0f(z,τ,m)
    Hux=[D1(Hu(m),x,t) for m in 1:M1]; Huxx=[D2(Hu(m),x,t) for m in 1:M1]
    Hut=[D1((z,τ)->Hc(z,τ)*ForwardDiff.derivative(b->u0f(z,b,m),τ),x,t) for m in 1:M1]
    ujhx=[D1((z,τ)->u0f(z,τ,m)*ForwardDiff.derivative(hbed,z),x,t) for m in 1:M1]
    ujHx=[D1((z,τ)->u0f(z,τ,m)*D1(Hc,z,τ),x,t) for m in 1:M1]
    etax = ForwardDiff.derivative(z->etaf_c(z,t), x)
    Lv(m) = (-ut[m]*hx, ut[m]*Hx, -Hut[m])
    Nv(k,j) = ( -u[j]*Huxx[k], D1((z,τ)->D1(Hu(k),z,τ)*u0f(z,τ,j),x,t),
                -u[k]*ujhx[j], u[k]*ujHx[j], -u[k]*Huxx[j],
                -Hux[j]*u[k]*hx/H, Hux[j]*u[k]*Hx/H, -Hux[j]*Hux[k]/H )
    acc = 0.0
    for j in 1:M1; acc += H*ut[j]*V.Mmat[i,j]; end
    for k in 1:M1, j in 1:M1
        acc += H*V.Mcal[i,k,j]*(u[k]*ux[j]) + Hux[k]*u[j]*V.Gcal[i,k,j]
    end
    acc += 9.81*H*etax*V.Phi[i]
    #  leading pressure: +∇( H²[ Σ L·P + Σ N·Pcal ] )   (moved to the LHS)
    LP(z,τ) = begin
        Hl=Hc(z,τ); Hxl=D1(Hc,z,τ); hxl=ForwardDiff.derivative(hbed,z)
        ul=[u0f(z,τ,m) for m in 1:M1]
        utl=[ForwardDiff.derivative(b->u0f(z,b,m),τ) for m in 1:M1]
        Hul(m)=(a,b)->Hc(a,b)*u0f(a,b,m)
        Huxl=[D1(Hul(m),z,τ) for m in 1:M1]; Huxxl=[D2(Hul(m),z,τ) for m in 1:M1]
        Hutl=[D1((a,b)->Hc(a,b)*ForwardDiff.derivative(e->u0f(a,e,m),b),z,τ) for m in 1:M1]
        ujhxl=[D1((a,b)->u0f(a,b,m)*ForwardDiff.derivative(hbed,a),z,τ) for m in 1:M1]
        ujHxl=[D1((a,b)->u0f(a,b,m)*D1(Hc,a,b),z,τ) for m in 1:M1]
        s1=zero(Hl)
        for j in 1:M1
            Lj=(-utl[j]*hxl, utl[j]*Hxl, -Hutl[j])
            for q in 1:3; s1 += Lj[q]*V.P[i,j,q]; end
        end
        for k in 1:M1, j in 1:M1
            Nkj=(-ul[j]*Huxxl[k], D1((a,b)->D1(Hul(k),a,b)*u0f(a,b,j),z,τ),
                 -ul[k]*ujhxl[j], ul[k]*ujHxl[j], -ul[k]*Huxxl[j],
                 -Huxl[j]*ul[k]*hxl/Hl, Huxl[j]*ul[k]*Hxl/Hl, -Huxl[j]*Huxl[k]/Hl)
            for q in 1:8; s1 += Nkj[q]*V.Pcal[i,k,j,q]; end
        end
        Hl^2*s1
    end
    acc += D1(LP, x, t)
    #  slope-pressure packages (moved to the LHS with a minus)
    for j in 1:M1
        Lj=Lv(j)
        for q in 1:3; acc -= H*(hx*Lj[q]*V.A[i,j,q] + Hx*Lj[q]*V.K[i,j,q]); end
    end
    for k in 1:M1, j in 1:M1
        Nkj=Nv(k,j)
        for q in 1:8; acc -= H*(hx*Nkj[q]*V.Acal[i,k,j,q] + Hx*Nkj[q]*V.Kcal[i,k,j,q]); end
    end
    acc
end

# ===========================================================================
#  STAGE 3 — our horizontal momentum residual vs Yang & Liu's PUBLISHED <2.30>
#
#  <2.30>:  sum_{n=0..4} ( R_{k-1,n} F1(k,n) + R_{k,n} F2(k,n) ) = 0
#  with R_{k,n} from supplementary section C (C.1)-(C.5) and F1,F2 from <2.28>,<2.29>.
#  This is the verification of the FULL NON-LINEAR model against the literature:
#  it uses their published coefficients, not the governing equation.
#
#  ⚠ TWO CONVENTIONS HAD TO BE RESOLVED, BOTH RECORDED IN README.md:
#  (a) their published p_k,n are the coefficients of -p_nh/(rho H), not +p_nh --
#      this follows from their own <3.4>, and section C's tails have the form of
#      the TRUE pressure polynomial, so section C must be fed with -p_k,n;
#  (b) the surface condition <2.5> already supplies the hydrostatic share of P_0
#      through P_1 = ... - g, so P_0 = C alone. Writing P_0 = g + C double-counts
#      gravity and shows up as a constant offset of exactly g*H_x in R_{k,0}.
# ===========================================================================
const GRAV = 9.81

"their w_{k,n}, <2.19>-<2.20> + (A.1)-(A.2)"
function wcoefs(c,x,t)
    M=length(c)-1
    H=Hc(x,t); Hx=D1(Hc,x,t); hx=ForwardDiff.derivative(hbed,x)
    u=[u0f(x,t,m) for m in 1:M+1]; ux=[D1((z,tau)->u0f(z,tau,m),x,t) for m in 1:M+1]
    w1=[(hx*u[k+1]-hx*u[k]+c[k+1]*H*ux[k]-c[k]*H*ux[k+1])/(c[k]-c[k+1]) for k in 1:M]
    w2=[(Hx*u[k]-Hx*u[k+1]+H*ux[k+1]-H*ux[k])/(2*(c[k]-c[k+1]))        for k in 1:M]
    w0=Vector{typeof(H)}(undef,M); w0[1]=-u[1]*hx
    for k in 2:M
        sm=-u[1]*hx
        for m in 1:k-1; sm+=c[m+1]*(w1[m]-w1[m+1])+c[m+1]^2*(w2[m]-w2[m+1]); end
        w0[k]=sm
    end
    w0,w1,w2
end
wn(c,x,t,k,n)=(W=wcoefs(c,x,t); n==0 ? W[1][k] : n==1 ? W[2][k] : W[3][k])

"their p_{k,n}, n=1..4, supplementary section B (B.1)-(B.4)"
function pcoefs_nh(c,x,t)
    M=length(c)-1
    H=Hc(x,t); Hx=D1(Hc,x,t); Ht=ForwardDiff.derivative(z->Hc(x,z),t)
    hx=ForwardDiff.derivative(hbed,x); u=[u0f(x,t,m) for m in 1:M+1]
    w0,w1,w2=wcoefs(c,x,t)
    w0x=[D1((z,tau)->wn(c,z,tau,k,0),x,t) for k in 1:M]
    w1x=[D1((z,tau)->wn(c,z,tau,k,1),x,t) for k in 1:M]
    w2x=[D1((z,tau)->wn(c,z,tau,k,2),x,t) for k in 1:M]
    w0t=[ForwardDiff.derivative(z->wn(c,x,z,k,0),t) for k in 1:M]
    w1t=[ForwardDiff.derivative(z->wn(c,x,z,k,1),t) for k in 1:M]
    w2t=[ForwardDiff.derivative(z->wn(c,x,z,k,2),t) for k in 1:M]
    P=Matrix{typeof(H)}(undef,M,4)
    for k in 1:M
        D=c[k]-c[k+1]
        P[k,1]= c[k]*hx*u[k+1]*w1[k]/(D*H) - c[k+1]*hx*u[k]*w1[k]/(D*H) +
                c[k]*u[k+1]*w0x[k]/D - c[k+1]*u[k]*w0x[k]/D + w0[k]*w1[k]/H + w0t[k]
        P[k,2]= hx*u[k]*w1[k]/(2D*H) - hx*u[k+1]*w1[k]/(2D*H) +
                c[k]*hx*u[k+1]*w2[k]/(D*H) - c[k+1]*hx*u[k]*w2[k]/(D*H) +
                c[k+1]*Hx*u[k]*w1[k]/(2D*H) - c[k]*Hx*u[k+1]*w1[k]/(2D*H) +
                u[k]*w0x[k]/(2D) + c[k]*u[k+1]*w1x[k]/(2D) -
                u[k+1]*w0x[k]/(2D) - c[k+1]*u[k]*w1x[k]/(2D) -
                Ht*w1[k]/(2H) + w1[k]^2/(2H) + w0[k]*w2[k]/H + w1t[k]/2
        P[k,3]= 2*hx*u[k]*w2[k]/(3D*H) - 2*hx*u[k+1]*w2[k]/(3D*H) +
                2*c[k+1]*Hx*u[k]*w2[k]/(3D*H) - 2*c[k]*Hx*u[k+1]*w2[k]/(3D*H) +
                Hx*u[k+1]*w1[k]/(3D*H) - Hx*u[k]*w1[k]/(3D*H) +
                u[k]*w1x[k]/(3D) + c[k]*u[k+1]*w2x[k]/(3D) -
                u[k+1]*w1x[k]/(3D) - c[k+1]*u[k]*w2x[k]/(3D) -
                2*Ht*w2[k]/(3H) + w1[k]*w2[k]/H + w2t[k]/3
        P[k,4]= Hx*u[k+1]*w2[k]/(2D*H) - Hx*u[k]*w2[k]/(2D*H) +
                u[k]*w2x[k]/(4D) - u[k+1]*w2x[k]/(4D) + w2[k]^2/(2H)
    end
    P
end

"full pressure polynomial P_0..P_4: sign convention (a) + surface condition (b)"
function Pfull(c,x,t)
    M=length(c)-1; p=pcoefs_nh(c,x,t)
    Pc=Matrix{eltype(p)}(undef,M,5)
    for k in 1:M
        Pc[k,2]=-p[k,1]-GRAV; Pc[k,3]=-p[k,2]; Pc[k,4]=-p[k,3]; Pc[k,5]=-p[k,4]
    end
    C=Vector{eltype(p)}(undef,M)
    C[M]=-sum(Pc[M,n+1] for n in 1:4)                       # <2.5>: p(sigma=1)=0
    for k in M-1:-1:1                                        # interface continuity, downward
        C[k]=C[k+1]+sum(Pc[k+1,n+1]*c[k+1]^n for n in 1:4)-sum(Pc[k,n+1]*c[k+1]^n for n in 1:4)
    end
    for k in 1:M; Pc[k,1]=C[k]; end
    Pc
end
Pn(c,x,t,k,n)=Pfull(c,x,t)[k,n+1]

"their R_{k,n}, supplementary section C (C.1)-(C.5)"
function Rcoefs(c,x,t)
    M=length(c)-1
    H=Hc(x,t); Hx=D1(Hc,x,t); hx=ForwardDiff.derivative(hbed,x)
    u=[u0f(x,t,m) for m in 1:M+1]; ux=[D1((z,tau)->u0f(z,tau,m),x,t) for m in 1:M+1]
    ut=[ForwardDiff.derivative(z->u0f(x,z,m),t) for m in 1:M+1]
    Ht=ForwardDiff.derivative(z->Hc(x,z),t)
    w0,w1,w2=wcoefs(c,x,t); Pc=Pfull(c,x,t)
    Px=[D1((z,tau)->Pn(c,z,tau,k,n),x,t) for k in 1:M, n in 0:4]
    R=Matrix{typeof(H)}(undef,M,5)
    for k in 1:M
        a=c[k]-c[k+1]; b=c[k+1]-c[k]
        R[k,1]=-c[k]*hx*u[k]*u[k+1]/(a*b*H)-c[k]*hx*u[k+1]^2/(b^2*H)-c[k+1]*hx*u[k]^2/(a^2*H)-
                c[k+1]*hx*u[k]*u[k+1]/(a*b*H)+u[k]*w0[k]/(a*H)+u[k+1]*w0[k]/(b*H)+
                c[k]^2*u[k+1]*ux[k+1]/b^2+c[k+1]*c[k]*u[k+1]*ux[k]/(a*b)+
                c[k+1]*c[k]*u[k]*ux[k+1]/(a*b)-c[k]*ut[k+1]/b+c[k+1]^2*u[k]*ux[k]/a^2-
                c[k+1]*ut[k]/a + hx*Pc[k,2]+Hx*Pc[k,1]+H*Px[k,1]
        R[k,2]= hx*u[k]^2/(a^2*H)+2*hx*u[k+1]*u[k]/(a*b*H)+hx*u[k+1]^2/(b^2*H)+
                c[k+1]*Hx*u[k]^2/(a^2*H)+c[k]*Hx*u[k+1]*u[k]/(a*b*H)+c[k+1]*Hx*u[k+1]*u[k]/(a*b*H)-
                Ht*u[k]/(a*H)+c[k]*Hx*u[k+1]^2/(b^2*H)-Ht*u[k+1]/(b*H)+
                u[k]*w1[k]/(a*H)+u[k+1]*w1[k]/(b*H)-2*c[k+1]*ux[k]*u[k]/a^2-
                c[k]*ux[k+1]*u[k]/(a*b)-c[k+1]*ux[k+1]*u[k]/(a*b)+ut[k]/a+ut[k+1]/b-
                c[k]*u[k+1]*ux[k]/(a*b)-c[k+1]*u[k+1]*ux[k]/(a*b)-2*c[k]*u[k+1]*ux[k+1]/b^2+
                2*hx*Pc[k,3]+H*Px[k,2]
        R[k,3]=-Hx*u[k]^2/(a^2*H)-2*Hx*u[k+1]*u[k]/(a*b*H)-Hx*u[k+1]^2/(b^2*H)+
                u[k]*w2[k]/(a*H)+u[k+1]*w2[k]/(b*H)+ux[k]*u[k]/a^2+ux[k+1]*u[k]/(a*b)+
                u[k+1]*ux[k]/(a*b)+u[k+1]*ux[k+1]/b^2+3*hx*Pc[k,4]-Hx*Pc[k,3]+H*Px[k,3]
        R[k,4]= 4*hx*Pc[k,5]-2*Hx*Pc[k,4]+H*Px[k,4]
        R[k,5]= H*Px[k,5]-3*Hx*Pc[k,5]
    end
    R
end

# <2.28>, <2.29>: the analytic weighting integrals
F1(c,k,n)=1/(c[k]-c[k-1])*((c[k]^(n+2)-c[k-1]^(n+2))/(n+2)-c[k-1]*(c[k]^(n+1)-c[k-1]^(n+1))/(n+1))
F2(c,k,n)=1/(c[k]-c[k+1])*((c[k+1]^(n+2)-c[k]^(n+2))/(n+2)-c[k+1]*(c[k+1]^(n+1)-c[k]^(n+1))/(n+1))

"<2.30>: the weighted momentum residual at each node"
weighted_theirs(c,x,t)=(M=length(c)-1; R=Rcoefs(c,x,t);
    [sum((k>1 ? R[k-1,n+1]*F1(c,k,n) : 0.0)+(k<M+1 ? R[k,n+1]*F2(c,k,n) : 0.0) for n in 0:4)
     for k in 1:M+1])

function stage3()
    println("STAGE 3 - our momentum residual vs Yang & Liu <2.30> (section C + <2.28>,<2.29>)\n")
    #  free check: F1+F2 at n=0 must equal the depth weight Phi_k
    c=[0.0,0.728,1.0]
    for k in 1:3
        f=(k>1 ? F1(c,k,0) : 0.0)+(k<3 ? F2(c,k,0) : 0.0)
        @printf("  F1+F2 vs Phi, node %d: %.12f  %.12f\n", k, f, Phi_w(c,k))
    end
    println()
    @printf("%-12s %-5s %-5s %-3s %16s %16s %10s\n","mesh","amp","x","k","ours","H x theirs","rel")
    println("-"^76)
    worst=0.0
    for cc in ([0.0,0.728,1.0],[0.0,0.4,0.8,1.0])
        CMESH[]=cc; V=assemble_vertical_tensors(length(cc)-1,1,cc)
        for e in (1.0,0.1), x in (3.7,12.0,19.4)
            AMP[]=e; th=weighted_theirs(cc,x,T0); H=Hc(x,T0)
            for i in 1:length(cc)
                a=Res_ours(V,cc,x,T0,i); b=H*th[i]
                rel=abs(a-b)/max(abs(b),1e-12); worst=max(worst,rel)
                @printf("M=%d nodes    %-5.2f %-5.1f %-3d %16.8e %16.8e %10.2e\n",
                        length(cc),e,x,i,a,b,rel)
            end
        end
    end
    AMP[]=1.0
    @printf("\nSTAGE 3 WORST RELATIVE DIFFERENCE: %.3e\n", worst)
end

run_stage2()
stage3()
