#!/bin/julia
include("constants.jl"); using constants
include("input.jl"); using input
using PyPlot
using Dierckx

dir = String[]

PyCall.PyDict(matplotlib["rcParams"])["font.size"] = "16"
PyCall.PyDict(matplotlib["rcParams"])["axes.titlesize"] = "24"
PyCall.PyDict(matplotlib["rcParams"])["axes.labelsize"] = "24"
PyCall.PyDict(matplotlib["rcParams"])["lines.linewidth"] = "3"
PyCall.PyDict(matplotlib["rcParams"])["lines.markersize"] = "8"
PyCall.PyDict(matplotlib["rcParams"])["legend.fontsize"] = "20"
PyCall.PyDict(matplotlib["rcParams"])["figure.subplot.left"] = "0.2"
PyCall.PyDict(matplotlib["rcParams"])["figure.subplot.bottom"] = "0.15"
PyCall.PyDict(matplotlib["rcParams"])["figure.subplot.wspace"] = "0.4"
PyCall.PyDict(matplotlib["rcParams"])["figure.subplot.hspace"] = "0.4"

function getdata(varname)
  global dir
  return readdlm(dir*"/"*varname*".dat")
end

#plt[:style][:use]("myfig")

read_input()

if length(ARGS) == 1
  dir = ARGS[1]
else
  dir = "test"
end

dir = "base"

if true
# Obtain all base data
  v= vec(getdata("vgrid"))
  r= vec(getdata("rgrid"))
  rgrid_gs2 = vec(getdata("rgrid_gs2"))
  rgrid_global = vec(getdata("rgrid_global"))

  Nrad = length(r)
  Nv = length(v)

  # 2D functions
  f0alpha = getdata("f0alpha")
  f0sd = getdata("f0sd")
  f0ash = getdata("f0ash")
  f0hot = getdata("f0hot")
  Drr = getdata("Drr")
  Drv = getdata("Drv")
  Dvr = getdata("Dvr")
  Dvv = getdata("Dvv")
  rflux = getdata("rflux")
  vflux = getdata("vflux")
  source = getdata("source")
  nus = getdata("nus")
  nupar = getdata("nupar")
  nuparv2 = getdata("nuparv2")
  totheatingv = getdata("totheatingv")
  iheatingv = getdata("iheatingv")
  eheatingv = getdata("eheatingv")
  sdiheatingv = getdata("sdiheatingv")
  sdeheatingv = getdata("sdeheatingv")
  sdpureiheatingv = getdata("sdpureiheatingv")
  sdpureeheatingv = getdata("sdpureeheatingv")
  dfdv = getdata("dfdv")
  dfdr = getdata("dfdr")

  # Functions of v only
  d3v = vec(getdata("d3v"))
  losF0 = vec(getdata("losF0"))
  losFsd = vec(getdata("losFsd"))
  
  # Functions of r
  reactionrate = vec(getdata("reactionrate"))
  nhot = vec(getdata("nhot"))
  nash = vec(getdata("nash"))
  nsd = vec(getdata("nsd"))
  nalpha = vec(getdata("nalpha"))
  Tash = vec(getdata("Tash"))
  Talpha = vec(getdata("Talpha"))
  pflux = vec(getdata("pflux"))
  hflux = vec(getdata("hflux"))
  flux0 = vec(getdata("flux0"))
  Dreff = vec(getdata("Dr_eff"))
  totheating = vec(getdata("totheating"))
  ionheating = vec(getdata("ionheating"))
  elheating = vec(getdata("elheating"))
  ashheating = vec(getdata("ashheating"))
  hotheating = vec(getdata("hotheating"))
  ihotheating = vec(getdata("ihotheating"))
  ehotheating = vec(getdata("ehotheating"))
  analheating = vec(getdata("analheating"))
  sdiheating = vec(getdata("sdiheating"))
  sdeheating = vec(getdata("sdeheating"))
  ashfit = vec(getdata("ashfit"))
  taus = vec(getdata("taus"))
  vcrit = vec(getdata("vcrit"))
  hfluxash = vec(getdata("hfluxash"))

  # Functions of global r
  Te_global = vec(getdata("Te"))
  Ti_global = vec(getdata("Ti"))
  ne_global = vec(getdata("ne"))*1.e20
  sourcetot = vec(getdata("sourcetot"))
  area_global = vec(getdata("area"))
  

  # Functions of GS2 rgrid
  chii = vec(getdata("chii"))
  phi2 = vec(getdata("phi2"))
  hflux_tot = vec(getdata("hflux_tot"))
  rhostar = vec(getdata("rhostar"))
end

function getgradient(y)
 temp = copy(y)
 temp[1] = (y[2] - y[1])/(r[2]-r[1])
 temp[end] = (y[end] - y[end-1])/(r[end]-r[end-1])
 temp[2:end-1] = (y[3:end]-y[1:end-2])./(r[3:end]-r[1:end-2])
 return temp
end

energy = 0.5*m_trace*v.^2

area_func = Spline1D(rgrid_global,area_global)
area= evaluate(area_func,r)

chii_func = Spline1D(rgrid_gs2,chii)
chii_rgrid = evaluate(chii_func,r)

Vprime_in = surface_area_in ./ grho_in
Vp_func = Spline1D(rgrid_in,Vprime_in)
Vprime = evaluate(Vp_func,r)

ne_func = Spline1D(rgrid_global,ne_global)
ne = evaluate(ne_func,r)

fsdpure = Array(Float64,(Nrad,Nv))
nsdpure = Array(Float64,Nrad)
sdpureiheating = Array(Float64,Nrad)
sdpureeheating = Array(Float64,Nrad)
Tsd = Array(Float64,Nrad)
for ir in 1:Nrad
  source_local = dot(d3v,vec(source[ir,:]))
  for iv in 1:Nv
    if (v[iv] <= valpha)
      fsdpure[ir,iv] = (0.25*source_local*taus[ir]/pi)/(vcrit[ir]^3+v[iv]^3)
    else
      fsdpure[ir,iv] = 0.0
    end
  end
  nsdpure[ir] = dot(d3v,vec(fsdpure[ir,:]))
  sdpureiheating[ir] = dot(d3v,vec(sdpureiheatingv[ir,:])./(4.0*pi*v.^2))
  sdpureeheating[ir] = dot(d3v,vec(sdpureeheatingv[ir,:])./(4.0*pi*v.^2))
  Tsd[ir] = dot(d3v.*energy,vec(f0sd[ir,:]))/nsd[ir]
end

#print(fsdpure[10,200])
#print(" ?= ")
#println(f0sd[10,200])

#if true

dir = "hirres"
r_hirres= vec(getdata("rgrid"))
totheating_hirres = vec(getdata("totheating"))
hotheating_hirres = vec(getdata("hotheating"))
ashheating_hirres = vec(getdata("ashheating"))
losF0_hirres = vec(getdata("losF0"))

dir = "hivres"
v_hivres= vec(getdata("vgrid"))
totheating_hivres = vec(getdata("totheating"))
hotheating_hivres = vec(getdata("hotheating"))
ashheating_hivres = vec(getdata("ashheating"))
losF0_hivres = vec(getdata("losF0"))


dir = "nedge1e18"
iheatingv_hin = getdata("iheatingv")
ionheating_hin = vec(getdata("ionheating"))

dir = "nedge1e18_hirres"
iheatingv_hirres = getdata("iheatingv")
ionheating_hirres = vec(getdata("ionheating"))

dir = "nedge1e18_hivres"
iheatingv_hivres = getdata("iheatingv")
ionheating_hivres = vec(getdata("ionheating"))

iheatingv_hivres_local = vec(iheatingv_hivres[1,:])
iheatingv_hirres_local = vec(iheatingv_hirres[1,:])
iheatingv_local = vec(iheatingv_hin[1,:])


iheatingv_hivres_local = vec(iheatingv_hivres[1,:])
iheatingv_hirres_local = vec(iheatingv_hirres[1,:])
iheatingv_local = vec(iheatingv_hin[1,:])

plt[:figure](figsize=(16,12))
plt[:tight_layout]
plt[:subplot](2,2,1)
plot(r/a,ionheating_hin*1.e-3,lw=1)
plot(r_hirres/a,ionheating_hirres*1.e-3,":")
plot(r/a,ionheating_hivres*1.e-3,"--")
ylabel(L"Ion heating rate $\left(\mathrm{kW}/\mathrm{m}^3\right)$")
xlabel(L"\psi / \psi_a")
#legend(("Std. res.",L"$N_\psi = 60$",L"$N_v = 800$"),loc="lower right")

plt[:subplot](2,2,2)
plot(v/valpha,iheatingv_local,lw=1)
plot(v/valpha,iheatingv_hirres_local,":")
plot(v_hivres/valpha,iheatingv_hivres_local,":")
ylabel("Integrand of ion heating rate")
xlabel(L"v / v_\alpha")
#legend(("Std. res.",L"$N_\psi = 60$",L"$N_v = 800$"),loc="lower right")


plt[:subplot](2,2,3)
semilogy(v/valpha,abs(losF0)/length(r),lw=1)
semilogy(v/valpha,abs(losF0_hirres)/length(r_hirres),":")
semilogy(v_hivres/valpha,abs(losF0_hivres)/length(r),"--")
ylabel(L"Line-of-sight averaged $F_0$ (A.U.)")
xlabel(L"$v / v_\alpha$")
legend(("Std. res.",L"$N_\psi = 60$",L"$N_v = 800$"),loc="upper right")

plt[:subplot](2,2,4)
plot(v/valpha,abs(losF0)/length(r),lw=1)
plot(v/valpha,abs(losF0_hirres)/length(r_hirres),":")
plot(v_hivres/valpha,abs(losF0_hivres)/length(r),"--")
#ylabel(L"Line-of-sight averaged $F_0$ (A.U.)")
xlabel(L"$v / v_\alpha$")
ylim(0.0,3e-5)
#legend(("Std. res.",L"$N_\psi = 60$",L"$N_v = 800$"),loc="upper right")
savefig("converged.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure]()
plt[:tight_layout]
plot(r/a,totheating,lw=1)
plot(r_hirres/a,totheating_hirres,":")
plot(r/a,totheating_hivres,"--")
ylabel(L"Total $\alpha$ heating rate $\left(\mathrm{W}/\mathrm{m}^3\right)$")
xlabel(L"\psi / \psi_a")
legend(("Std. res.",L"$N_\psi = 60$",L"$N_v = 800$"),loc="upper right")
savefig("conv_heatv.png",bbox_inches="tight")
cla()
clf()
close()


Tifunc = Spline1D(rgrid_global,Ti_global)
Ti_gs2 = evaluate(Tifunc,rgrid_gs2)
Ti = evaluate(Tifunc,r)
Tefunc = Spline1D(rgrid_global,Te_global)
Te_gs2 = evaluate(Tefunc,rgrid_gs2)
Te= evaluate(Tefunc,r)
cs = sqrt(Te_gs2/mref)
vti_gs2 = sqrt(2.0*Ti_gs2/mref)
vti= sqrt(2.0*evaluate(Tifunc,r)/mref)
rhos = a*rhostar.*cs./vti_gs2
chiGB_gs2 = rhos.^2.*cs/a
chifunc = Spline1D(rgrid_gs2,chiGB_gs2,k=1)
chiGB = evaluate(chifunc,r)
cA_0p5 = vti_gs2[1]/sqrt(0.00682)
cA_0p65 = 0.5*(vti_gs2[2]+vti_gs2[3])/sqrt(0.0047)

println(chiGB_gs2)



ir_sample = int(Nrad/2)
plt[:figure](figsize=(16,6))
plt[:subplot](1,2,1)
plt[:tight_layout]
semilogy(v/valpha,abs(vec(f0alpha[ir_sample,:])),"-k")
semilogy(v/valpha,abs(vec(f0sd[ir_sample,:])),"--r")
ylim(1.e-6,0.1)
xlim(0.0,1.1)
vlines(vcrit[ir_sample]/valpha,1e-6,2e-6)
text(vcrit[ir_sample]/valpha,2e-6,L"$v_c$",fontsize=32)
#vlines(cA_0p65/valpha,1e-6,2e-6)
#text(cA_0p65/valpha,2e-6,L"$c_A$",fontsize=32)
ylabel(L"$F_{0\alpha}$ at $\psi = 0.65 \psi_a$ (A.U.)")
xlabel(L"$v / v_\alpha$")
legend((L"Numerical $F_0$","Analytic slowing-down + ash"),fontsize=14)

plt[:subplot](1,2,2)

plot(v/valpha,abs(vec(f0alpha[ir_sample,:])),"-k")
plot(v/valpha,abs(vec(f0sd[ir_sample,:])),"--r")
ylim(0.0,4e-5)
xlim(0.0,1.1)
xlabel(L"$v / v_\alpha$")

savefig("f0sd.png",bbox_inches="tight")
cla()
clf()
close()

ir_sample = 1
plt[:figure]()
plt[:tight_layout]
#ir_sample = int(length(r)/2)
ir_sample = 1
semilogy(v/valpha,abs(vec(f0alpha[ir_sample,:])),"-k")
semilogy(v/valpha,abs(vec(f0sd[ir_sample,:])),"--r")
ylim(1.e-6,0.1)
xlim(0.0,1.1)
vlines(vcrit[ir_sample]/valpha,1e-6,2e-6)
text(vcrit[ir_sample]/valpha,2e-6,L"$v_c$",fontsize=32)
#vlines(cA_0p5/valpha,1e-6,2e-6)
#text(cA_0p5/valpha,2e-6,L"$c_A$",fontsize=32)
ylabel(L"$F_{0\alpha}$ at $\psi = 0.5 \psi_a$ (A.U.)")
xlabel(L"$v / v_\alpha$")
legend((L"Numerical $F_0$","Analytic slowing-down + ash"),fontsize=20)
savefig("f0sd-midrad.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure]()
plt[:tight_layout]
ir_sample = int(length(r)/3)
#ir_sample = 1
semilogy(v/valpha,abs(vec(f0alpha[ir_sample,:])),"-k")
semilogy(v/valpha,abs(vec(f0sd[ir_sample,:])),"--r")
ylim(1.e-6,0.1)
xlim(0.0,1.1)
vlines(vcrit[ir_sample]/valpha,1e-6,2e-6)
text(vcrit[ir_sample]/valpha,2e-6,L"$v_c$",fontsize=32)
#vlines(cA_0p5/valpha,1e-6,2e-6)
#text(cA_0p5/valpha,2e-6,L"$c_A$",fontsize=32)
ylabel(L"$F_{0\alpha}$ at $\psi = 0.6 \psi_a$ (A.U.)")
xlabel(L"$v / v_\alpha$")
legend((L"Numerical $F_0$","Analytic slowing-down + ash"),fontsize=20)
savefig("f0sd-0p6.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure](figsize=(6,2))
plt[:tight_layout]
ir_sample = int(length(r)/3)
#ir_sample = 1
semilogy(v/valpha,abs(vec(f0alpha[ir_sample,:]))*1e4,"-k",lw=2)
ir_sample = 2*ir_sample
semilogy(v/valpha,abs(vec(f0alpha[ir_sample,:]))*1e4,"-c",lw=2)
#vlines(cA_0p5/valpha,1e-6,2e-6)
#text(cA_0p5/valpha,2e-6,L"$c_A$",fontsize=32)

#ir_sample = 1
ir_sample = int(length(r)/3)
semilogy(v/valpha,abs(vec(f0sd[ir_sample,:]))*1e4,"--k",lw=2)
ir_sample = 2*ir_sample
semilogy(v/valpha,abs(vec(f0sd[ir_sample,:]))*1e4,"--c",lw=2)
ylim(1.e-2,10.0)
xlim(0.0,1.1)
#vlines(vcrit[ir_sample]/valpha,1e-6,2e-6)
#text(vcrit[ir_sample]/valpha,2e-6,L"$v_c$",fontsize=32)
ylabel(L"$F_{0\alpha}$ (A.U.)",fontsize=18)
xlabel(L"$v / v_\alpha$",fontsize=18)
#legend((L"Numerical $F_0$","Analytic slowing-down + ash"),fontsize=20)
legend((L"$r / a = $ 0.6",L"$r / a = $ 0.7"),fontsize=14)
savefig("f0sd-multi.png",bbox_inches="tight")
cla()
clf()
close()



Tifunc = Spline1D(rgrid_global,Ti_global)
Ti_gs2 = evaluate(Tifunc,rgrid_gs2)
Ti = evaluate(Tifunc,r)
Tefunc = Spline1D(rgrid_global,Te_global)
Te_gs2 = evaluate(Tefunc,rgrid_gs2)
Te= evaluate(Tefunc,r)
cs = sqrt(Te_gs2/mref)
vti_gs2 = sqrt(2.0*Ti_gs2/mref)
vti= sqrt(2.0*evaluate(Tifunc,r)/mref)
rhos = a*rhostar.*cs./vti_gs2
chiGB_gs2 = rhos.^2.*cs/a
chifunc = Spline1D(rgrid_gs2,chiGB_gs2,k=1)
chiGB = evaluate(chifunc,r)
cA_0p5 = vti_gs2[1]/sqrt(0.00682)
cA_0p65 = 0.5*(vti_gs2[2]+vti_gs2[3])/sqrt(0.0047)

Dvv_phys = Array(Float64,(Nrad,Nv))
for ir in 1:Nrad
  Drr[ir,:] = Drr[ir,:]/chiGB[ir]
  Drv[ir,:] = (a/valpha)*Drv[ir,:]/chiGB[ir]
  Dvr[ir,:] = (a/valpha)*Dvr[ir,:]/chiGB[ir]
  Dvv_phys[ir,:] = Dvv[ir,:]
  Dvv[ir,:] = (a/valpha)^2*Dvv[ir,:]/chiGB[ir]
end

plt[:figure](figsize=(16,12))
plt[:tight_layout]
plt[:subplot](2,2,1)
    pcolormesh(r/a,v/valpha,Drr')
    title(L"$D_{\psi \psi} / \chi_\mathrm{GB}$"*"\n",fontsize=24)
    ylim(0.0,1.0)
    xlim(0.5,0.8)
    ylabel(L"v / v_\alpha",fontsize=32)
    colorbar()
plt[:tight_layout]
plt[:subplot](2,2,2)
    pcolormesh(r/a,v/valpha,Drv',vmin=-0.5,vmax=0.0)
    title(L"$D_{\psi v} a / v_\alpha \chi_\mathrm{GB}$"*"\n",fontsize=24)
    ylim(0.0,1.0)
    xlim(0.5,0.8)
    colorbar()
plt[:tight_layout]
plt[:subplot](2,2,3)
    pcolormesh(r/a,v/valpha,Dvr',vmin=-0.5,vmax=0.0)
    title(L"$D_{v \psi} a / v_\alpha \chi_\mathrm{GB}$"*"\n",fontsize=24)
    ylabel(L"v / v_\alpha",fontsize=32)
    xlabel(L"\psi / \psi_a",fontsize=32)
    ylim(0.0,1.0)
    xlim(0.5,0.8)
    colorbar()
plt[:tight_layout]
plt[:subplot](2,2,4)
    pcolormesh(r/a,v/valpha,Dvv',vmin=0.0,vmax=0.05)
    title(L"$D_{v v} a^2 / v_\alpha^2 \chi_\mathrm{GB}$"*"\n",fontsize=24)
    ylim(0.0,1.0)
    xlim(0.5,0.8)
    xlabel(L"\psi / \psi_a",fontsize=32)
    colorbar()
plt[:tight_layout]
savefig("diffcontour.png",bbox_inches="tight")
cla()
clf()
close()

ir = 20
plt[:figure]()
plt[:tight_layout]
ylabel(L"Normalized diffusion coefficient at (Abs) $\left(\mathrm{m}^2/s\right)$",fontsize=16)
xlabel(L"$v / v_\alpha$")
semilogy(v/valpha,abs(vec(Drr[ir,:])),"-k")
semilogy(v/valpha,abs(vec(Drv[ir,:])),"--g")
semilogy(v/valpha,abs(vec(Dvr[ir,:])),"-.r")
semilogy(v/valpha,abs(vec(Dvv[ir,:])),":b")
legend((L"$D_{\psi \psi}$",L"$D_{\psi v}$",L"$D_{v \psi}$",L"$D_{v v}$"),fontsize=16)
savefig("difflog.png",bbox_inches="tight")
cla()
clf()
close()


ir_sample=1
dir = "constdiff0p1"
f0_D0 = vec(getdata("f0alpha")[ir_sample,:])
dir = "constdiff1"
f0_D1 = vec(getdata("f0alpha")[ir_sample,:])
dir = "constdiff3"
f0_D3 = vec(getdata("f0alpha")[ir_sample,:])
dir = "constdiff6"
f0_D6 = vec(getdata("f0alpha")[ir_sample,:])
dir = "constdiff15"
f0_D15 = vec(getdata("f0alpha")[ir_sample,:])
dir = "constdiff30"
f0_D30 = vec(getdata("f0alpha")[ir_sample,:])

plt[:figure]()
plt[:tight_layout]
#semilogy(energy/(1.e6*el),vec(abs(f0sd[ir_sample,:])))
plot(energy/(1.e6*el),abs(f0_D0),"-")
plot(energy/(1.e6*el),abs(f0_D1),"--")
plot(energy/(1.e6*el),abs(f0_D3),"-.")
plot(energy/(1.e6*el),abs(f0_D6),":")
plot(energy/(1.e6*el),abs(f0_D15),"--")
#ylim(1.e-8,1.e-5)
ylim(0.0,1.e-5)
xlim(0.0,3.7)
ylabel(L"$F_0$ at mid-radius (A.U.)")
xlabel(L"$E \left(\mathrm{MeV}\right)$ ")
#legend((L"$D=0.1 \mathrm{m}^2/s$",L"$D=1\mathrm{m}^2/s$",L"$D=3\mathrm{m}^2/s$",L"$D=6\mathrm{m}^2/s$",L"$D=15\mathrm{m}^2/s$",L"$D=30\mathrm{m}^2/s$"),fontsize=12,loc="lower left")
legend((L"$D=0.1 \mathrm{m}^2/\mathrm{s}$",L"$D=1 \mathrm{m}^2/\mathrm{s}$",L"$D=3\mathrm{m}^2/\mathrm{s}$",L"$D=6\mathrm{m}^2/\mathrm{s}$",L"$D=15\mathrm{m}^2/\mathrm{s}$",L"$D=30\mathrm{m}^2/\mathrm{s}$"),fontsize=12,loc="upper right")
savefig("sigmarcomp.png",bbox_inches="tight")
cla()
clf()
close()

dir = "weakturb"
rflux_weak = getdata("rflux")
vflux_weak = getdata("vflux")
totheat_weak = getdata("totheating")
ionheat_weak = getdata("ionheating")
nalpha_weak = getdata("nalpha")
elheat_weak = getdata("elheating")
dir = "strongturb"
losF0_strong = getdata("losF0")
rflux_strong = getdata("rflux")
vflux_strong = getdata("vflux")
totheat_strong = getdata("totheating")
ionheat_strong = getdata("ionheating")
elheat_strong = getdata("elheating")
nalpha_strong = getdata("nalpha")

edgespectrum = 2.0*pi*m_trace*vec(rflux[end,:]).*v.^4
edgespectrum_weak = 2.0*pi*m_trace*vec(rflux_weak[end,:]).*v.^4
edgespectrum_strong = 2.0*pi*m_trace*vec(rflux_strong[end,:]).*v.^4

plt[:figure](figsize=(6,2))
plt[:tight_layout]
plot(v/valpha,edgespectrum,"-k")
plot(v/valpha,edgespectrum_strong,"--r")
plot(v/valpha,edgespectrum_weak,":b")
legend(("Nominal amplitude",L"$\chi_i \times$ 5",L"$\chi_i \times$ 0.2",L"$\chi_i \times$ 0.01"),loc="upper right",fontsize=12)
xlabel(L"$v / v_\alpha$",fontsize=18)
ylabel(L"$2 \pi m_\alpha v^2 \Gamma_r \, \left(\mathrm{J} / \mathrm{m}^3 \right)$",fontsize=14)
xlim(0.0,1.0)
ylim(-0.004,0.015)
savefig("edgespectrum.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure](figsize=(6,2))
plt[:tight_layout]
semilogy(v/valpha,abs(edgespectrum),"-k")
semilogy(v/valpha,abs(edgespectrum_strong),"--r")
semilogy(v/valpha,abs(edgespectrum_weak),":b")
#legend(("Nominal amplitude",L"$\chi_i \times$ 5",L"$\chi_i \times$ 0.2",L"$\chi_i \times$ 0.01"),loc="upper right",fontsize=12)
xlabel(L"$v / v_\alpha$",fontsize=18)
ylabel(L"$|2 \pi m_\alpha v^2 \Gamma_r \, \left(\mathrm{J} / \mathrm{m}^3 \right)|$",fontsize=14)
savefig("edgespectrum_log.png",bbox_inches="tight")
cla()
clf()
close()


#if false
#for ir in 1:Nrad
#  for iv in 1:Nv
#    rflux[ir,iv] = rflux[ir,iv]*Vprime[ir]
#    rflux_weak[ir,iv] = rflux_weak[ir,iv]*Vprime[ir]
#    rflux_strong[ir,iv] = rflux_strong[ir,iv]*Vprime[ir]
#    vflux[ir,iv] = vflux[ir,iv]*v[iv]^2
#    vflux_weak[ir,iv] = vflux_weak[ir,iv]*v[iv]^2
#    vflux_strong[ir,iv] = vflux_strong[ir,iv]*v[iv]^2
#  end
#end
#end

plt[:tight_layout]
streamplot(r/a,v/valpha,rflux',vflux'*a/valpha)
#streamplot(r/a,v/valpha,rflux',vflux'*a/valpha)
#streamplot(r/a,v/valpha,rflux',vflux',color=sqrt(rflux.^2+(vflux).^2)')
#streamplot(r,v,rflux',vflux')
xlabel(L"$\psi / \psi_a$",fontsize=36)
ylabel(L"$v / v_\alpha$",fontsize=36)
xlim(0.5,0.8)
ylim(0.0,1.05)
savefig("basestream.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure](figsize=(25,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
streamplot(r/a,v/valpha,rflux_weak',vflux_weak'*a/valpha)
xlabel(L"$\psi / \psi_a$",fontsize=36)
#streamplot(r/a,v/valpha,rflux_weak',vflux_weak',color=sqrt(rflux_weak.^2+(vflux_weak).^2)')
title(L"$\chi_{i,\mathrm{nominal}} / 5$"*"\n",fontsize=24)
xlim(0.5,0.8)
ylim(0.0,1.0)

plt[:subplot](1,3,2)
streamplot(r/a,v/valpha,rflux',vflux'*a/valpha)
#streamplot(r/a,v/valpha,rflux',vflux',color=sqrt(rflux.^2+(vflux).^2)')
title(L"$\chi_{i,\mathrm{nominal}}$"*"\n",fontsize=24)
xlim(0.5,0.8)
ylim(0.0,1.0)
xlabel(L"$\psi / \psi_a$",fontsize=36)

plt[:subplot](1,3,3)
streamplot(r/a,v/valpha,rflux_strong',vflux_strong'*(a/valpha))
#streamplot(r/a,v/valpha,rflux_strong',vflux_strong',color=sqrt(rflux_strong.^2+(vflux_strong).^2)')
title(L"$\chi_{i,\mathrm{nominal}} \times 5$"*"\n",fontsize=24)
xlim(0.5,0.8)
ylim(0.0,1.0)
xlabel(L"$\psi / \psi_a$",fontsize=36)

savefig("streamcomp.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure]()
plt[:tight_layout]
plot(r/a,totheat_weak,":g")
plot(r/a,totheating,"-k")
plot(r/a,totheat_strong,"--r")
legend((L"$\chi_{i,\mathrm{nominal}}/5$",L"$\chi_{i,\mathrm{nominal}}$",L"$\chi_{i,\mathrm{nominal}}\times 5$"),fontsize=36)
xlabel(L"$\psi / \psi_a$",fontsize=48)
ylabel(L"Total $\alpha$ heating $\left( \mathrm{W}/\mathrm{m}^3 \right)$",fontsize=48)
savefig("strengthcomp.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure]()
plt[:tight_layout]
plot(rgrid_gs2/a,phi2,"-ok")
ylabel(L"$\sum_\mathbf{k} | e \phi_\mathbf{k} / T_i |^2$",fontsize=36)
xlim(0.45,0.85)
xlabel(L"$\psi / \psi_a$",fontsize=32)
savefig("phi2.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure]()
ax1 = plt[:axes]()
#ax1.plot(r/a,chii,"-b")
#ax1.set_xlabel(L"$\psi / \psi_a$")
#ax1.set_ylabel(L"$\chi_i / \chi_\mathrm{GB}$")
plot(rgrid_gs2/a,chii./chiGB_gs2,"-ok")
xlabel(L"$\psi / \psi_a$",fontsize=32)
xlim(0.45,0.85)
ylim(0.0,11.0)
ax1[:yaxis][:label][:set_color]("black")
ylabel(L"$\chi_i / \chi_\mathrm{GB}$",fontsize=32)
for t1 in ax1[:get_yticklabels]()
  t1[:set_color]("k")
end

ax2 = plt[:twinx]()
#ax2=[:axes]()
#ax2.plot(r/a,chii,"Db")
#ax2.set_ylabel(L"$\chi_i$ $\mathrm{m}^2/s$$")
plot(rgrid_gs2/a,chii,"--Dg")
ylabel(L"$\chi_i$ $\left(\mathrm{m}^2/\mathrm{s}\right)$",fontsize=24)
ax2[:yaxis][:label][:set_color]("green")
for t1 in ax2[:get_yticklabels]()
  t1[:set_color]("g")
end
xlim(0.45,0.85)
ylim(0.0,11.0)
savefig("chii.png",bbox_inches="tight")
cla()
clf()
close()

dir="biasin"
chii_in = vec(getdata("chii"))
totheat_in = vec(getdata("totheating"))
ionheat_in = vec(getdata("ionheating"))
elheat_in = vec(getdata("elheating"))
nalpha_in = vec(getdata("nalpha"))
dir="biasout"
chii_out = vec(getdata("chii"))
totheat_out = vec(getdata("totheating"))
ionheat_out = vec(getdata("ionheating"))
elheat_out = vec(getdata("elheating"))
nalpha_out = vec(getdata("nalpha"))

plt[:figure]()
plot(rgrid_gs2/a,chii_in./chiGB_gs2,"--vb")
plot(rgrid_gs2/a,chii./chiGB_gs2,"-ok")
plot(rgrid_gs2/a,chii_out./chiGB_gs2,"--^r")
ylabel(L"$\chi_i / \chi_\mathrm{GB}$")
xlabel(L"$\psi / \psi_a$")
legend(("Biased in","Nominal", "Biased out"),loc="upper left")
plt[:tight_layout]
savefig("biaschii.png",bbox_inches="tight")
cla()
clf()
close()


plt[:figure]()
plt[:tight_layout]
plot(r/a,totheat_in,":g")
plot(r/a,totheating,"-k")
plot(r/a,totheat_out,"--r")
legend(("Biased in",L"Nominal $\chi_i$","Biased out"))
xlabel(L"$\psi / \psi_a$",fontsize=36)
ylabel(L"Total $\alpha$ heating $\left( \mathrm{W}/\mathrm{m}^3 \right)$")
savefig("biascomp.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure]()
plt[:tight_layout]
plot(r/a,totheating,"-k")
plot(r/a,ionheating,"--b")
plot(r/a,elheating,"--g")
plot(r/a,ashheating,":",color="0.5")
legend((L"Total $\alpha$ heating","Ion heating","Electron heating","Heating from ash"))
ylabel(L"Heating rate $\left( \mathrm{W}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a$",fontsize=36)
savefig("speciesheatcomp.png",bbox_inches="tight")
cla()
clf()
close()

dir="nedge1e18"
totheat_1e18 = vec(getdata("totheating"))
ionheat_1e18 = vec(getdata("ionheating"))
elheat_1e18 = vec(getdata("elheating"))
ashheat_1e18 = vec(getdata("ashheating"))
analheat_1e18 = vec(getdata("analheating"))
heatv_1e18 = getdata("totheatingv")
iheatv_1e18 = getdata("iheatingv")
eheatv_1e18 = getdata("eheatingv")
rflux_1e18 = getdata("rflux")
hfluxash_1e18 = vec(getdata("hfluxash"))
dir="nedge1e19"
totheat_1e19 = vec(getdata("totheating"))
ionheat_1e19 = vec(getdata("ionheating"))
elheat_1e19 = vec(getdata("elheating"))
ashheat_1e19 = vec(getdata("ashheating"))
analheat_1e19 = vec(getdata("analheating"))

plt[:figure]()
plt[:tight_layout]
plot(r/a,totheating,"-k")
plot(r/a,totheat_1e18,"-g")
plot(r/a,totheat_1e19,"-r")
plot(r/a,ashheating,"--k")
plot(r/a,ashheat_1e18,"--g")
plot(r/a,ashheat_1e19,"--r")
legend((L"$n_{\alpha,\mathrm{edge}} = 10^{17}/\mathrm{m}^3$",L"$n_{\alpha,\mathrm{edge}} = 10^{18}/\mathrm{m}^3$",L"$n_{\alpha,\mathrm{edge}} = 10^{19}/\mathrm{m}^3$"),fontsize=12)
ylabel(L"Heating rate $\left( \mathrm{W}/\mathrm{m}^3 \right)$")
xlabel(L"$\psi / \psi_a$",fontsize=36)
savefig("nedgecomp.png",bbox_inches="tight")
cla()
clf()
close()

dir="eject"
nalpha_eject = vec(getdata("nalpha"))
totheat_eject = vec(getdata("totheating"))
rflux_eject = getdata("rflux")
vflux_eject = getdata("vflux")


plt[:figure]()
plt[:tight_layout]
plot(r/a,nalpha,"k")
plot(r/a,nalpha_eject,"r")
xlabel(L"$\psi / \psi_a$",fontsize=36)
ylabel(L"$n_\alpha$ $\left(1/\mathrm{m}^3\right)$")
legend(("Nominal case","Ejected core"))
savefig("ejectcomp_nalpha.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure]()
plt[:tight_layout]
plot(r/a,totheating,"k")
plot(r/a,totheat_eject,"r")
xlabel(L"$\psi / \psi_a$",fontsize=36)
ylabel(L"Total $\alpha$ heating rate $\left( \mathrm{W}/\mathrm{m}^3 \right)$")
legend(("Nominal case","Ejected core"))
savefig("ejectcomp_heat.png",bbox_inches="tight")
cla()
clf()
close()



dir="circular"
nalpha_circ = vec(getdata("nalpha"))
totheat_circ = vec(getdata("totheating"))
ionheat_circ = vec(getdata("totheating"))

plt[:figure]()
plt[:tight_layout]
plot(r/a,totheating,"k")
plot(r/a,totheat_circ,"r")
xlabel(L"$\psi / \psi_a$",fontsize=36)
ylabel(L"$\alpha$-ion heating rate $\left( \mathrm{W}/\mathrm{m}^3 \right)$")
legend(("Nominal case","Circular geometry"),fontsize=16)
savefig("circcomp.png",bbox_inches="tight")
cla()
clf()
close()

dir="noeflux"
nalpha_noeflux = vec(getdata("nalpha"))
totheat_noeflux = vec(getdata("totheating"))
ionheat_noeflux = vec(getdata("ionheating"))
dir="diffonly"
nalpha_diffonly = vec(getdata("nalpha"))
totheat_diffonly = vec(getdata("totheating"))
ionheat_diffonly = vec(getdata("ionheating"))

plt[:figure]()
plt[:tight_layout]
plot(r/a,totheating,"k")
plot(r/a,totheat_noeflux,"--g")
plot(r/a,totheat_diffonly,":r")
xlabel(L"$\psi / \psi_a$",fontsize=36)
ylabel(L"Total $\alpha$ heating rate $\left( \mathrm{W}/\mathrm{m}^3 \right)$")
legend(("Nominal case",L"$\Gamma_v=0$",L"$D_{\psi \psi}$ only"),fontsize=16)
savefig("diffmodel_totheat.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure]()
plt[:tight_layout]
plot(r/a,ionheating,"k")
plot(r/a,ionheat_noeflux,"--g")
plot(r/a,ionheat_diffonly,":r")
xlabel(L"$\psi / \psi_a$",fontsize=36)
ylabel(L"$\alpha$-ion heating rate $\left( \mathrm{W}/\mathrm{m}^3 \right)$")
legend(("Nominal case",L"No $\Gamma_v$",L"$D_{\psi \psi}$ only"),fontsize=16)
savefig("diffmodel_ionheat.png",bbox_inches="tight")
cla()
clf()
close()


plt[:figure]()
plt[:tight_layout]
plot(r/a,nalpha,"k")
plot(r/a,nalpha_noeflux,"--g")
plot(r/a,nalpha_diffonly,":r")
xlabel(L"$\psi / \psi_a$",fontsize=36)
ylabel(L"$n_\alpha$ $\left(1/\mathrm{m}^3\right)$")
legend(("Nominal case",L"No $\Gamma_v$",L"$D_{\psi \psi}$ only"),fontsize=16)
savefig("diffmodel_nalpha.png",bbox_inches="tight")
cla()
clf()
close()

ir=20
plt[:figure]()
plt[:tight_layout]
semilogy(v/valpha,abs(vec(nuparv2[ir,:])),"m")
semilogy(v/valpha,abs(vec(Dvv_phys[ir,:])),"c")
xlabel(L"$v / v_\alpha$",fontsize=36)
ylabel("Energy diffusion (A.U.)")
legend((L"$\frac{1}{2}\nu_{\|\|} v^2 $",L"$D_{vv}$"),fontsize=16)
savefig("nuparVSDvv.png",bbox_inches="tight")
cla()
clf()
close()


plt[:figure]()
plt[:tight_layout]
semilogy(energy/(1.e6*el),abs(losF0_strong)/length(r),"-k")
semilogy(energy/(1.e6*el),abs(losFsd)/length(r),"--r")
ylim(1.e-6,0.1)
xlim(0.0,3.7)
ylabel(L"Line-of-sight averaged $F_0$ (A.U.)")
xlabel(L"$E \left(\mathrm{MeV}\right)$ ")
legend((L"Numerical $F_0$","Analytic slowing-down + ash"),fontsize=20)
savefig("LOSstrong.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure]()
plt[:tight_layout]
#semilogy(energy/(1.e6*el),abs(losF0)/length(r),"-k")
#semilogy(energy/(1.e6*el),abs(losFsd)/length(r),"--r")
semilogy(v/valpha,abs(losF0)/length(r),"-k")
semilogy(v/valpha,abs(losFsd)/length(r),"--r")
ylim(1.e-6,0.1)
xlim(0.0,1.1)
ylabel(L"Line-of-sight averaged $F_0$ (A.U.)")
xlabel(L"$v / v_\alpha$",fontsize=32)
legend((L"Numerical $F_0$","Analytic slowing-down + ash"),fontsize=20)
savefig("LOS.png",bbox_inches="tight")
cla()
clf()
close()


plt[:figure](figsize=(12,6))
plot(r/a,nalpha./ne,"-b")
plot(r/a,nhot./ne,"-.r")
plot(r/a,nsdpure./ne,"--g")
plot(r/a,nash./ne,":k")
ylabel(L"$n / n_e$",fontsize=36)
xlabel(L"$\psi / \psi_a$",fontsize=36)
legend((L"$n_\alpha$",L"$n_\mathrm{hot}$",L"$n_\mathrm{sd}$",L"$n_\mathrm{ash}$"),fontsize=24)
savefig("nalphaVSr.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure](figsize=(16,6))
plt[:subplot](1,2,1)
plt[:tight_layout]
plot(r/a,totheating*1.e-3,"-b")
plot(r/a,(sdiheating+sdeheating)*1.e-3,"--r")
ylabel(L"$\alpha$ heating rate $\left(\mathrm{kW} / \mathrm{m}^3\right)$",fontsize=24)
xlabel(L"$\psi / \psi_a$",fontsize=36)
#legend((L"$n_\alpha$",L"$n_\mathrm{hot}$",L"$n_\mathrm{sd}$",L"$n_\mathrm{ash}$"),fontsize=24)
#legend((L"$n_\alpha$",L"$n_\mathrm{sd}$"),fontsize=24)
plt[:subplot](1,2,2)
plt[:tight_layout]
plot(r/a,area.*totheating*1.e-3,"-b")
plot(r/a,area.*(sdiheating+sdeheating)*1.e-3,"--r")
ylabel("Area-weighted \n"*L"$\alpha$ heating rate $\left(\mathrm{kW} / \mathrm{m} \right)$",fontsize=24)
xlabel(L"$\psi / \psi_a$",fontsize=36)
legend((L"$F_{0\alpha}$ from T3CORE","Analytic slowing-down"),fontsize=18)
savefig("alphaheat.png",bbox_inches="tight")
cla()
clf()
close()

ir = 1
plt[:figure](figsize=(16,6))
plt[:subplot](1,2,1)

plot(r/a,ionheating[:]*1.e-3,"b")
plot(r/a,elheating[:]*1.e-3,"g")
plot(r/a,sdiheating[:]*1.e-3,"b--")
plot(r/a,sdeheating[:]*1.e-3,"g--")
ylabel(L"Heating rate by species $\left( \mathrm{kW}/ \mathrm{m}^3 \right)$",fontsize=24)
xlabel(L"$\psi / \psi_a$",fontsize=36)
legend(("Ions","Electrons"),fontsize=24)

plt[:subplot](1,2,2)
plot(v/valpha,vec(iheatingv[ir,:]),"b")
plot(v/valpha,vec(eheatingv[ir,:]),"g")
#plot(v/valpha,vec(sdpureiheatingv[ir,:]),"--b")
#plot(v/valpha,vec(sdpureeheatingv[ir,:]),"--g")
plot(v/valpha,vec(sdiheatingv[ir,:]),"--b")
plot(v/valpha,vec(sdeheatingv[ir,:]),"--g")
text(vcrit[ir]/valpha,-.028,L"$v_c$",fontsize=32)
vlines(vcrit[ir_sample]/valpha,-.03,-.028)
ylabel("Integrand of heating rate",fontsize=24)
xlabel(L"$v / v_\alpha$",fontsize=36)
ylim(-.03,.03)
#xlim(-.02,.02)
savefig("heatrv.png",bbox_inches="tight")
cla()
clf()
close()


ir = 15
plt[:figure](figsize=(16,6))
plt[:subplot](1,2,1)

plot(r/a,ionheat_1e18[:],"b")
plot(r/a,elheat_1e18[:],"g")
ylabel(L"Heating rate by species $\left( \mathrm{W}/ \mathrm{m}^3 \right)$",fontsize=24)
xlabel(L"$\psi / \psi_a$",fontsize=36)
legend(("Ions","Electrons"),fontsize=24)

plt[:subplot](1,2,2)
plot(v/valpha,vec(iheatv_1e18[ir,:]),"b")
plot(v/valpha,vec(eheatv_1e18[ir,:]),"g")
text(vcrit[ir]/valpha,-.028,L"$v_c$",fontsize=32)
vlines(vcrit[ir_sample]/valpha,-.03,-.028)
ylabel("Integrand of heating rate",fontsize=24)
xlabel(L"$v / v_\alpha$",fontsize=36)
ylim(-.03,.03)
#xlim(-.02,.02)
savefig("heatrv_1e18.png",bbox_inches="tight")
cla()
clf()
close()

ir = 15
integrand = 4.0*pi*(v.^2).*energy.*vec(rflux[ir,:])
integrand_1e18 = 4.0*pi*(v.^2).*energy.*vec(rflux_1e18[ir,:])
plt[:figure]()
plot(v/valpha,integrand,"-k")
plot(v/valpha,integrand_1e18,"--r")
ylabel("Integrand of heat flux",fontsize=24)
xlabel(L"$v / v_\alpha$",fontsize=36)
legend(("Nominal",L"$n_\mathrm{edge} = 10^{18}$"))
ylim(-.03,.03)
#xlim(-.02,.02)
savefig("hfluxintegrand.png",bbox_inches="tight")
cla()
clf()
close()





plt[:figure](figsize=(16,6))
plt[:subplot](1,2,1)

plot(r/a,ionheating[:],"b")
plot(r/a,elheating[:],"g")
plot(r/a,sdpureiheating[:],"b--")
plot(r/a,sdpureeheating[:],"g--")
ylabel(L"Heating rate by species $\left( \mathrm{W}/ \mathrm{m}^3 \right)$",fontsize=24)
xlabel(L"$\psi / \psi_a$",fontsize=36)
legend(("Ions","Electrons"),fontsize=24)


plt[:subplot](1,2,2)
plot(v/valpha,vec(iheatingv[ir,:]),"b")
plot(v/valpha,vec(eheatingv[ir,:]),"g")
plot(v/valpha,vec(sdpureiheatingv[ir,:]),"--b")
plot(v/valpha,vec(sdpureeheatingv[ir,:]),"--g")
text(vcrit[ir]/valpha,-.028,L"$v_c$",fontsize=32)
vlines(vcrit[ir_sample]/valpha,-.03,-.028)
ylabel("Integrand of heating rate",fontsize=24)
xlabel(L"$\psi / \psi_a$",fontsize=36)
ylim(-.03,.03)
#xlim(-.02,.02)
savefig("heatrvpure.png",bbox_inches="tight")
cla()
clf()
close()

dir = "nedge1e16"
fit_1e16 = vec(getdata("ashfit"))
nash_1e16 = vec(getdata("nash"))
Tash_1e16 = vec(getdata("Tash"))
totheat_1e16 = vec(getdata("totheating"))
eheat_1e16 = vec(getdata("elheating"))
iheat_1e16 = vec(getdata("ionheating"))
ashheat_1e16 = vec(getdata("ashheating"))
hotheat_1e16 = vec(getdata("hotheating"))
dir = "nedge5e16"
fit_5e16 = vec(getdata("ashfit"))
nash_5e16 = vec(getdata("nash"))
Tash_5e16 = vec(getdata("Tash"))
totheat_5e16 = vec(getdata("totheating"))
eheat_5e16 = vec(getdata("elheating"))
iheat_5e16 = vec(getdata("ionheating"))
ashheat_5e16 = vec(getdata("ashheating"))
hotheat_5e16 = vec(getdata("hotheating"))
hflux_5e16 = vec(getdata("hflux"))
rflux_5e16 = getdata("rflux")
dir = "nedge8e16"
fit_8e16 = vec(getdata("ashfit"))
nash_8e16 = vec(getdata("nash"))
Tash_8e16 = vec(getdata("Tash"))
totheat_8e16 = vec(getdata("totheating"))
eheat_8e16 = vec(getdata("elheating"))
iheat_8e16 = vec(getdata("ionheating"))
ashheat_8e16 = vec(getdata("ashheating"))
hotheat_8e16 = vec(getdata("hotheating"))
fit_1e17 = copy(ashfit)
dir = "nedge2e17"
fit_2e17 = vec(getdata("ashfit"))
nash_2e17 = vec(getdata("nash"))
Tash_2e17 = vec(getdata("Tash"))
totheat_2e17 = vec(getdata("totheating"))
eheat_2e17 = vec(getdata("elheating"))
iheat_2e17 = vec(getdata("ionheating"))
ashheat_2e17 = vec(getdata("ashheating"))
hotheat_2e17 = vec(getdata("hotheating"))
dir = "nedge3e17"
fit_3e17 = vec(getdata("ashfit"))
nash_3e17 = vec(getdata("nash"))
Tash_3e17 = vec(getdata("Tash"))
totheat_3e17 = vec(getdata("totheating"))
eheat_3e17 = vec(getdata("elheating"))
iheat_3e17 = vec(getdata("ionheating"))
ashheat_3e17 = vec(getdata("ashheating"))
hotheat_3e17 = vec(getdata("hotheating"))
dir = "nedge5e17"
fit_5e17 = vec(getdata("ashfit"))
nash_5e17 = vec(getdata("nash"))
Tash_5e17 = vec(getdata("Tash"))
totheat_5e17 = vec(getdata("totheating"))
eheat_5e17 = vec(getdata("elheating"))
iheat_5e17 = vec(getdata("ionheating"))
ashheat_5e17 = vec(getdata("ashheating"))
hotheat_5e17 = vec(getdata("hotheating"))
hflux_5e17 = vec(getdata("hflux"))
dir = "nedge1e18"
fit_1e18 = vec(getdata("ashfit"))
nash_1e18 = vec(getdata("nash"))
Tash_1e18 = vec(getdata("Tash"))
totheat_1e18 = vec(getdata("totheating"))
eheat_1e18 = vec(getdata("elheating"))
iheat_1e18 = vec(getdata("ionheating"))
ashheat_1e18 = vec(getdata("ashheating"))
hotheat_1e18 = vec(getdata("hotheating"))
hflux_1e18 = vec(getdata("hflux"))
f0alpha_1e18 = getdata("f0alpha")
f0ash_1e18 = getdata("f0ash")
#dir = "nedge2e18"

ir = 1
plt[:figure](figsize=(16,6))
plt[:subplot](1,2,1)
plot(v/vti[ir],vec(f0alpha_1e18[ir,:]),"+b")
plot(v/vti[ir],vec(f0ash_1e18[ir,:]),"-k")
xlim(0.0,1.0)
xlabel(L"$v / v_{ti}$")
ylabel(L"Distribution function $F_{0}$ (A.U.)")
legend((L"$F_{0\alpha}$ from T3CORE",L"$F_\mathrm{ash}$"))

plt[:subplot](1,2,2)
#plot(r/a,fit_1e16)
plot(r/a,fit_5e16,"Dr")
#plot(r/a,fit_8e16)
plot(r/a,fit_1e17,"om")
#plot(r/a,fit_2e17)
plot(r/a,fit_5e17,"^g")
plot(r/a,fit_1e18,"vb")
xlabel(L"$\psi / \psi_a $")
ylabel("RMS error of fit")
#legend((L"n_\mathrm{edge} = 10^{16}/\mathrm{m}^3",L"n_\mathrm{edge} = 5 \times 10^{16}/\mathrm{m}^3",L"n_\mathrm{edge} = 10^{17}/\mathrm{m}^3",L"n_\mathrm{edge} = 5^{17}/\mathrm{m}^3",L"n_\mathrm{edge} = 10^{18}/\mathrm{m}^3"),fontsize=12)
legend((L"n_\mathrm{edge} = 5 \times 10^{16}/\mathrm{m}^3",L"n_\mathrm{edge} = 10^{17}/\mathrm{m}^3",L"n_\mathrm{edge} = 5 \times 10^{17}/\mathrm{m}^3",L"n_\mathrm{edge} = 10^{18}/\mathrm{m}^3"),fontsize=16,loc="lower right")
savefig("goodfit.png",bbox_inches="tight")
cla()
clf()
close()


plt[:figure](figsize=(16,6))
plt[:subplot](1,2,1)
#plot(r/a,nash_1e16)
plot(r/a,nash_5e16*1.e-18,"--r")
plot(r/a,nash*1.e-18,"-.m")
#plot(r/a,nash_2e17,":")
plot(r/a,nash_5e17*1.e-18,":g")
plot(r/a,nash_1e18*1.e-18,"-.b")
plot(r/a,ne*1.e-20,"-y")
text(0.55,0.95,L"$n_e/100$",color="y",fontsize=16)
xlabel(L"$\psi / \psi_a $")
ylabel(L"$n_\mathrm{ash} \, \left( 10^{18} / \mathrm{m}^{3}\right)$")
#legend((L"n_\mathrm{edge} = 5 \times 10^{16}/\mathrm{m}^3",L"n_\mathrm{edge} = 10^{17}/\mathrm{m}^3",L"n_\mathrm{edge} = 5 \times 10^{17}/\mathrm{m}^3",L"n_\mathrm{edge} = 10^{18}/\mathrm{m}^3"),fontsize=16)

plt[:subplot](1,2,2)
plot(r/a,Tash_5e16/(1000.0*el),"--r")
plot(r/a,Tash/(1000.0*el),"-.m")
#plot(r/a,Tash_2e17/(1000.0*el),":")
plot(r/a,Tash_5e17/(1000.0*el),":g")
plot(r/a,Tash_1e18/(1000.0*el),"-.b")
plot(r/a,Ti/(1000.0*el),"-k")
plot(r/a,Te/(1000.0*el),"-y")
xlabel(L"$\psi / \psi_a $")
ylabel(L"$T_\mathrm{ash} \left( \mathrm{keV} \right)$")
legend((L"n_\mathrm{edge} = 5 \times 10^{16}/\mathrm{m}^3",L"n_\mathrm{edge} = 10^{17}/\mathrm{m}^3",L"n_\mathrm{edge} = 5 \times 10^{17}/\mathrm{m}^3",L"n_\mathrm{edge} = 10^{18}/\mathrm{m}^3"),fontsize=16,loc="upper right")
text(0.55,14.75,L"$T_e$",color="y",fontsize=16)
text(0.6,11.25,L"$T_i$",color="k",fontsize=16)
savefig("ashprofiles.png",bbox_inches="tight")
cla()
clf()
close()

plot(r/a,Tash_5e16/(1000.0*el),"--r")
plot(r/a,Tash/(1000.0*el),"-.m")
plot(r/a,Tash_5e17/(1000.0*el),":g")
plot(r/a,Tash_1e18/(1000.0*el),"-.b")
plot(r/a,Ti/(1000.0*el),"-k")
plot(r/a,Te/(1000.0*el),"--k")
xlabel(L"$\psi / \psi_a $")
ylabel(L"$T_\mathrm{ash} \left( \mathrm{keV} \right)$")
legend((L"n_{\alpha\mathrm{edge}} = 5 \times 10^{16}/\mathrm{m}^3",L"n_{\alpha\mathrm{edge}} = 10^{17}/\mathrm{m}^3",L"n_{\alpha\mathrm{edge}} = 5 \times 10^{17}/\mathrm{m}^3",L"n_{\alpha\mathrm{edge}} = 10^{18}/\mathrm{m}^3"),fontsize=16,loc="upper right")
ylim(6.0,16.0)
text(0.55,14.75,L"$T_e$",color="k",fontsize=24)
text(0.6,11.25,L"$T_i$",color="k",fontsize=24)
savefig("ashtemp.png",bbox_inches="tight")
cla()
clf()
close()


plt[:figure](figsize=(16,12))
#plt[:subplot](2,2,1)
#plt[:tight_layout]
#plot(r/a,totheat_5e16,"--r")
#plot(r/a,totalheating,"-.m")
#plot(r/a,totheat_5e17,":g")
#plot(r/a,totheat_1e18,"-.b")

plt[:subplot](2,2,1)
plt[:tight_layout]
title("Ash")
plot(r/a,ashheat_5e16*1.e-3,"--r")
plot(r/a,ashheating*1.e-3,"-.m")
plot(r/a,ashheat_5e17*1.e-3,":g")
plot(r/a,ashheat_1e18*1.e-3,"-.b")
xlim(0.5,0.8)
ylim(-120.0,160.0)
ylabel(L"$\alpha$ heating rate $\left( \mathrm{kW}/ \mathrm{m}^3\right)$")


plt[:subplot](2,2,2)
plt[:tight_layout]
title(L"Fast $\alpha$ particles")
plot(r/a,hotheat_5e16*1.e-3,"--r")
plot(r/a,hotheating*1.e-3,"-.m")
plot(r/a,hotheat_5e17*1.e-3,":g")
plot(r/a,hotheat_1e18*1.e-3,"-.b")
xlim(0.5,0.8)
ylim(-120.0,160.0)

plt[:subplot](2,2,3)
plt[:tight_layout]
title("Ions")
plot(r/a,iheat_5e16*1.e-3,"--r")
plot(r/a,ionheating*1.e-3,"-.m")
plot(r/a,iheat_5e17*1.e-3,":g")
plot(r/a,iheat_1e18*1.e-3,"-.b")
xlim(0.5,0.8)
ylim(-120.0,160.0)
ylabel(L"$\alpha$ heating rate $\left( \mathrm{kW}/ \mathrm{m}^3 \right)$")
xlabel(L"$\psi / \psi_a $")

plt[:subplot](2,2,4)
plt[:tight_layout]
title("Electrons")
plot(r/a,eheat_5e16*1.e-3,"--r")
plot(r/a,elheating*1.e-3,"-.m")
plot(r/a,eheat_5e17*1.e-3,":g")
plot(r/a,eheat_1e18*1.e-3,"-.b")
xlim(0.5,0.8)
ylim(-120.0,160.0)
xlabel(L"$\psi / \psi_a $")
legend((L"n_\mathrm{edge} = 5 \times 10^{16}/\mathrm{m}^3",L"n_\mathrm{edge} = 10^{17}/\mathrm{m}^3",L"n_\mathrm{edge} = 5 \times 10^{17}/\mathrm{m}^3",L"n_\mathrm{edge} = 10^{18}/\mathrm{m}^3"),fontsize=16,loc="lower right")
savefig("allheating.png",bbox_inches="tight")
cla()
clf()
close()

gradTash = Array(Float64,Nrad)
gradTash[1] = (Tash[2] - Tash[1])/(r[2]-r[1])
gradTash[end] = (Tash[end] - Tash[end-1])/(r[end]-r[end-1])
gradTash[2:end-1] = (Tash[3:end]-Tash[1:end-2])./(r[3:end]-r[1:end-2])
chiash = -hfluxash./(gradTash.*nash)

gradTalpha = Array(Float64,Nrad)
gradTalpha[1] = (Talpha[2] - Talpha[1])/(r[2]-r[1])
gradTalpha[end] = (Talpha[end] - Talpha[end-1])/(r[end]-r[end-1])
gradTalpha[2:end-1] = (Talpha[3:end]-Talpha[1:end-2])./(r[3:end]-r[1:end-2])
chialpha = -hflux./(gradTalpha.*nalpha)

plt[:figure]()
plt[:tight_layout]
plot(r/a,chialpha./chii_rgrid,"-b")
plot(r/a,chiash./chii_rgrid,"--",color="0.5")
xlabel(L"$\psi / \psi_a $")
ylabel(L"$\chi_\alpha / \chi_i$")
xlim(0.5,0.8)
legend((L"Total $\alpha$","Ash only"))
savefig("chialpha.png",bbox_inches="tight")
cla()
clf()
close()

p0 = ne_global[1] * Te_global[1]

#####
# Slowing-down comparison
#####

plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,(sdiheating+sdeheating)*1.e-3,"--r")
legend(("With turbulence","Analytic slowing-down + ash"),fontsize=16)
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,getgradient(nsd.*Tsd)*a/p0,"--r")
#legend(("With turbulence","Analytic slowing-down"),fontsize=16)
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux*1.e-3,"-k")
#plot(r/a,hfluxash*1.e-3,"--r")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

savefig("sdcomp.png",bbox_inches="tight")
cla()
clf()
close()

#####
# Turbulence strength
#####

dir = "strongturb"
hflux_strong = vec(getdata("hflux"))
heat_strong = vec(getdata("totheating"))
T_strong = vec(getdata("Talpha"))
n_strong = vec(getdata("nalpha"))

dir = "weakturb"
hflux_weak = vec(getdata("hflux"))
heat_weak = vec(getdata("totheating"))
T_weak = vec(getdata("Talpha"))
n_weak = vec(getdata("nalpha"))

dir = "veryweakturb"
hflux_vweak = vec(getdata("hflux"))
heat_vweak = vec(getdata("totheating"))
T_vweak = vec(getdata("Talpha"))
n_vweak = vec(getdata("nalpha"))


plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,heat_strong*1.e-3,"--r")
plot(r/a,heat_weak*1.e-3,":b")
plot(r/a,heat_vweak*1.e-3,"-.c")
ylim(1.0,250.0)
#title("Heating",fontsize=36)
legend(("Nominal amplitude",L"$\chi_i \times$ 5",L"$\chi_i \times$ 0.2",L"$\chi_i \times$ 0.01"),loc="upper right",fontsize=18)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=36)
xlabel(L"$r / a $",fontsize=36)

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,-getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,-getgradient(n_strong.*T_strong)*a/p0,"--r")
plot(r/a,-getgradient(n_weak.*T_weak)*a/p0,"-.b")
plot(r/a,-getgradient(n_vweak.*T_vweak)*a/p0,"-.c")
ylim(0.1,0.99)
#title(L"$\nabla p_\alpha$",fontsize=36)
ylabel(L"$ -\left(\partial p_\alpha / \partial r \right) \, \times a / p_{e0}$",fontsize=36)
xlabel(L"$r / a $",fontsize=36)

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux*1.e-3,"-k")
plot(r/a,hflux_strong*1e-3,"--r")
plot(r/a,hflux_weak*1e-3,"-.b")
plot(r/a,hflux_vweak*1e-3,"-.c")
ylim(0.1,60.0)
#title("Heat flux",fontsize=36)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=36)
xlabel(L"$r / a $",fontsize=36)

savefig("strengthcomp.png",bbox_inches="tight")
cla()
clf()
close()

plt[:tight_layout]
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,heat_strong*1.e-3,"--r")
plot(r/a,heat_weak*1.e-3,":b")
plot(r/a,heat_vweak*1.e-3,"-.c")
ylim(1.0,250.0)
xlim(0.5,0.8)
#title("Heating",fontsize=36)
legend(("Nominal amplitude",L"$\chi_i \times$ 5",L"$\chi_i \times$ 0.2",L"$\chi_i \times$ 0.01"),loc="upper right",fontsize=18)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=36)
xlabel(L"$r / a $",fontsize=36)
savefig("strengthcomp-heating.png",bbox_inches="tight")
cla()
clf()
close()

plt[:tight_layout]
plot(r/a,-getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,-getgradient(n_strong.*T_strong)*a/p0,"--r")
plot(r/a,-getgradient(n_weak.*T_weak)*a/p0,"-.b")
plot(r/a,-getgradient(n_vweak.*T_vweak)*a/p0,"-.c")
ylim(0.1,0.99)
xlim(0.5,0.8)
#title(L"$\nabla p_\alpha$",fontsize=36)
ylabel(L"$ -\left(\partial p_\alpha / \partial r \right) \, \times a / p_{e0}$",fontsize=36)
xlabel(L"$r / a $",fontsize=36)
savefig("strengthcomp-pressure.png",bbox_inches="tight")
cla()
clf()
close()




plt[:tight_layout]
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,heat_strong*1.e-3,"--r")
plot(r/a,heat_weak*1.e-3,":b")
legend(("Nominal amplitude",L"$\times 5$",L"$\times 0.2$",L"$\times 0.01$"),fontsize=16)
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")
savefig("heatingstrength.png",bbox_inches="tight")
cla()
clf()
close()




plt[:tight_layout]
plot(r/a,getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,getgradient(n_strong.*T_strong)*a/p0,"--r")
plot(r/a,getgradient(n_weak.*T_weak)*a/p0,":b")
#plot(r/a,getgradient(n_vweak.*T_vweak)*a/p0,"-.c")
title("Pressure gradient",fontsize=36)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")
savefig("pressuregradientstrength.png",bbox_inches="tight")
cla()
clf()
close()



plt[:tight_layout]
plot(r/a,nalpha.*Talpha/p0,"-k")
plot(r/a,n_strong.*T_strong/p0,"--r")
plot(r/a,n_weak.*T_weak/p0,":b")
ylabel(L"$n_\alpha T_\alpha \, / n_{e0} T_{e0} $",fontsize=16)
xlabel(L"$\psi / \psi_a $")
title("Pressure profile",fontsize=36)
legend(("Nominal amplitude",L"$\chi_i \times 5$",L"$\chi_i / 5$"),fontsize=16)
savefig("pressureprofilestrength.png",bbox_inches="tight")
cla()
clf()
close()

plt[:tight_layout]
plot(r/a,hflux*1.e-3,"-k")
plot(r/a,hflux_strong*1e-3,"--r")
plot(r/a,hflux_weak*1e-3,":b")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")
legend(("Nominal amplitude",L"$\chi_i \times 5$",L"$\chi_i / 5$"),fontsize=16)
savefig("hfluxstrength.png",bbox_inches="tight")
cla()
clf()
close()

#####
# Bias
#####

dir = "biasout"
hflux_out = vec(getdata("hflux"))
heat_out = vec(getdata("totheating"))
T_out = vec(getdata("Talpha"))
n_out = vec(getdata("nalpha"))

dir = "biasin"
hflux_in = vec(getdata("hflux"))
heat_in = vec(getdata("totheating"))
T_in = vec(getdata("Talpha"))
n_in = vec(getdata("nalpha"))

plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,heat_out*1.e-3,"--r")
plot(r/a,heat_in*1.e-3,"-.b")
legend(("Nominal amplitude","Biased outward","Biased inward"),fontsize=16)
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,getgradient(n_out.*T_out)*a/p0,"--r")
plot(r/a,getgradient(n_in.*T_in)*a/p0,"-.b")
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux*1.e-3,"-k")
plot(r/a,hflux_out*1.e-3,"--r")
plot(r/a,hflux_in*1.e-3,"-.b")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

savefig("biascomp.png",bbox_inches="tight")
cla()
clf()
close()

#####
# edge density
#####

dir = "nedge1e18"
hflux_1e18 = vec(getdata("hflux"))
heat_1e18 = vec(getdata("totheating"))
iheat_1e18 = vec(getdata("ionheating"))
T_1e18 = vec(getdata("Talpha"))
n_1e18 = vec(getdata("nalpha"))

dir = "nedge5e17"
hflux_5e17 = vec(getdata("hflux"))
heat_5e17 = vec(getdata("totheating"))
T_5e17 = vec(getdata("Talpha"))
n_5e17 = vec(getdata("nalpha"))

dir = "nedge5e16"
hflux_5e16 = vec(getdata("hflux"))
heat_5e16 = vec(getdata("totheating"))
T_5e16 = vec(getdata("Talpha"))
n_5e16 = vec(getdata("nalpha"))


plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,heat_5e16*1.e-3,"-.b")
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,heat_5e17*1.e-3,"-.m")
plot(r/a,heat_1e18*1.e-3,"--r")
legend((L"n_\mathrm{edge} = 5 \times 10^{16}/\mathrm{m}^3",L"n_\mathrm{edge} = 10^{17}/\mathrm{m}^3",L"n_\mathrm{edge} = 5 \times 10^{17}/\mathrm{m}^3",L"n_\mathrm{edge} = 10^{18}/\mathrm{m}^3"),fontsize=16)
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(n_5e16.*T_5e16)*a/p0,":b")
plot(r/a,getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,getgradient(n_5e17.*T_5e17)*a/p0,"-.m")
plot(r/a,getgradient(n_1e18.*T_1e18)*a/p0,"--r")
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux_5e16*1.e-3,":b")
plot(r/a,hflux*1.e-3,"-k")
plot(r/a,hflux_5e17*1.e-3,"-.m")
plot(r/a,hflux_1e18*1.e-3,"--r")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

savefig("nedgecomp.png",bbox_inches="tight")
cla()
clf()
close()

#####
# diffusion model
#####

dir = "diffonly"
hflux_diff = vec(getdata("hflux"))
heat_diff = vec(getdata("totheating"))
T_diff = vec(getdata("Talpha"))
n_diff = vec(getdata("nalpha"))

dir = "noeflux"
hflux_noe = vec(getdata("hflux"))
heat_noe = vec(getdata("totheating"))
T_noe = vec(getdata("Talpha"))
n_noe = vec(getdata("nalpha"))

dir = "nedge1e18_diffonly"
hflux_diff1e18 = vec(getdata("hflux"))
heat_diff1e18 = vec(getdata("totheating"))
T_diff1e18 = vec(getdata("Talpha"))
n_diff1e18 = vec(getdata("nalpha"))



plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,heat_noe*1.e-3,"--r")
plot(r/a,heat_diff*1.e-3,":b")
legend(("Nominal",L"No $\Gamma_v$",L"$D_{\psi \psi}$ only"),fontsize=16)
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,getgradient(n_noe.*T_noe)*a/p0,"--r")
plot(r/a,getgradient(n_diff.*T_diff)*a/p0,":b")
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux*1.e-3,"-k")
plot(r/a,hflux_noe*1.e-3,"--r")
plot(r/a,hflux_diff*1.e-3,":b")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

savefig("diffmodelcomp.png",bbox_inches="tight")
cla()
clf()
close()



plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,heat_1e18*1.e-3,"--r")
plot(r/a,heat_diff1e18*1.e-3,":b")
legend(("Nominal",L"$n_\mathrm{edge} = 10^{18}/ \mathrm{m}^3$",L"$D_{\psi \psi}$ only"),fontsize=16)
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,getgradient(n_1e18.*T_1e18)*a/p0,"--r")
plot(r/a,getgradient(n_diff1e18.*T_diff1e18)*a/p0,":b")
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux*1.e-3,"-k")
plot(r/a,hflux_1e18*1.e-3,"--r")
plot(r/a,hflux_diff1e18*1.e-3,":b")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

savefig("hindiffmodelcomp.png",bbox_inches="tight")
cla()
clf()
close()

#####
# Ejection mode
#####

dir = "eject"
hflux_eject = vec(getdata("hflux"))
heat_eject = vec(getdata("totheating"))
T_eject = vec(getdata("Talpha"))
n_eject = vec(getdata("nalpha"))
dir = "eject_strong"
hflux_eject_strong = vec(getdata("hflux"))
heat_eject_strong = vec(getdata("totheating"))
T_eject_strong = vec(getdata("Talpha"))
n_eject_strong = vec(getdata("nalpha"))
dir = "eject_weak"
hflux_eject_weak = vec(getdata("hflux"))
heat_eject_weak = vec(getdata("totheating"))
T_eject_weak = vec(getdata("Talpha"))
n_eject_weak = vec(getdata("nalpha"))

plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,heat_eject_strong*1.e-3,"--r")
plot(r/a,heat_eject*1.e-3,"-.g")
plot(r/a,heat_eject_weak*1.e-3,":b")
legend(("Nominal",L"Direct ejection, $\chi_i \times 5$",L"Direct ejection, $\chi_{i,\mathrm{nominal}}$",L"Direct ejection, $\chi_i / 5$"),fontsize=16)
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,getgradient(n_eject_strong.*T_eject_strong)*a/p0,"--r")
plot(r/a,getgradient(n_eject.*T_eject)*a/p0,"-.g")
plot(r/a,getgradient(n_eject_weak.*T_eject_weak)*a/p0,":b")
ylim(-60.0,0.0)
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux*1.e-3,"-k")
plot(r/a,hflux_eject_strong*1.e-3,"--r")
plot(r/a,hflux_eject*1.e-3,"-.g")
plot(r/a,hflux_eject_weak*1.e-3,":b")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

savefig("ejectcomp.png",bbox_inches="tight")
cla()
clf()
close()

plt[:tight_layout]
plot(r/a,hflux*1.e-3,"-k")
#plot(r/a,hflux_eject_strong*1.e-3,"--r")
plot(r/a,hflux_eject*1.e-3,":b")
legend(("Nominal","Core ejection"))
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

savefig("hfluxeject.png",bbox_inches="tight")
cla()
clf()
close()


#####
# Analytic model 
#####

dir="esanal"
heat_esanal = vec(getdata("totheating"))
hflux_esanal = vec(getdata("hflux"))
n_esanal = vec(getdata("nalpha"))
T_esanal = vec(getdata("Talpha"))

plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,heat_esanal*1.e-3,"--r")
legend(("Nominal","Analytic model"),fontsize=16)
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,getgradient(n_esanal.*T_esanal)*a/p0,"--r")
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux*1.e-3,"-k")
plot(r/a,hflux_esanal*1.e-3,"--r")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

savefig("esanalcomp.png",bbox_inches="tight")
cla()
clf()
close()

#####
# Geometry
#####

dir="circular"
heat_circ= vec(getdata("totheating"))
hflux_circ = vec(getdata("hflux"))
n_circ = vec(getdata("nalpha"))
T_circ = vec(getdata("Talpha"))

plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,heat_circ*1.e-3,"--r")
legend(("Nominal","Circular geometry"),fontsize=16)
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,getgradient(n_circ.*T_circ)*a/p0,"--r")
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux*1.e-3,"-k")
plot(r/a,hflux_circ*1.e-3,"--r")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

savefig("circcomp.png",bbox_inches="tight")
cla()
clf()
close()


#####
# EM model (crap)
####

dir="emanal"
heat_emanal = vec(getdata("totheating"))
hflux_emanal = vec(getdata("hflux"))
n_emanal = vec(getdata("nalpha"))
T_emanal = vec(getdata("Talpha"))

plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,heat_esanal*1.e-3,"-.b")
plot(r/a,heat_emanal*1.e-3,"--r")
legend(("Nominal","ES analytic model","EM analytic model"),fontsize=16)
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,getgradient(n_esanal.*T_esanal)*a/p0,"-.b")
plot(r/a,getgradient(n_emanal.*T_emanal)*a/p0,"--r")
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux*1.e-3,"-k")
plot(r/a,hflux_esanal*1.e-3,"-.b")
plot(r/a,hflux_emanal*1.e-3,"--r")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

savefig("emanalcomp.png",bbox_inches="tight")
cla()
clf()
close()

#####
# Source
#####

rgrid_plot = linspace(rgrid_global[1],rgrid_global[end],1000)
source_func = Spline1D(rgrid_global,sourcetot)
sourcearea_func = Spline1D(rgrid_global,sourcetot.*area_global)
source_plot = evaluate(source_func,rgrid_plot)
sourcearea_plot = evaluate(sourcearea_func,rgrid_plot)

plt[:figure](figsize=(16,6))
plt[:subplot](1,2,1)
plt[:tight_layout]
plot(rgrid_plot/a,source_plot,"-g")
ylabel(L"$\alpha$ source $ \sigma_\alpha \left( 1/ \mathrm{m}^3 \mathrm{s} \right)  $")
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,2,2)
plt[:tight_layout]
plot(rgrid_plot/a,sourcearea_plot,"-g")
ylabel(L"Area-weighted source $\left( 1/ \mathrm{m}-\mathrm{s} \right)  $")
xlabel(L"$\psi / \psi_a $")
axvspan(0.5,0.8,facecolor="0.2",alpha=0.2)

savefig("sourcer.png",bbox_inches="tight")
cla()
clf()
close()

ir=10
println(size(source))
println(size(v))
plt[:figure](figsize=(16,6))
plt[:subplot](1,2,1)
plt[:tight_layout]
plot(v/valpha,vec(source[ir,:]) ,"-g")
ylabel(L"$\alpha source $")
xlabel(L"$v / v_\alpha $")

plt[:subplot](1,2,2)
plt[:tight_layout]
plot(rgrid_plot/a,sourcearea_plot,"-g")
ylabel(L"Area-weighted source $\left( 1/ \mathrm{m}-\mathrm{s} \right)  $")
xlabel(L"$\psi / \psi_a $")
axvspan(0.5,0.8,facecolor="0.2",alpha=0.2)

savefig("sourcerv.png",bbox_inches="tight")
cla()
clf()
close()



tau_eq = (1.8e-19 * sqrt(4.0 * 2.0) * 4.0 * 0.5 *1.0e14* 15.0)^(-1) * ((2.0*Tash + 4.0*Ti).^1.5)/el^1.5
tau_eq += (1.8e-19 * sqrt(4.0 * 3.0) * 4.0 * 0.5 *1.0e14* 15.0)^(-1) * ((3.0*Tash + 4.0*Ti).^1.5)/el^1.5
tau_eq += (1.8e-19 * sqrt(4.0 * (me/mp)) * 4.0 * 1.0 *1.0e14* 15.0)^(-1) * (((me/mp)*Tash + 4.0*Te).^1.5)/el^1.5

#print("tau_eq = "); println(tau_eq)

vref = sqrt(2.0*Ti[end]/(2.0*mp))
tau_conf = sum(Vprime_in.*(Te_global+Ti_global))*1.0e20*(rgrid_global[2]-rgrid_global[1])
tau_conf = tau_conf / (rhostar[end]^2*1.0e20*vref*Ti[end]*hflux_tot[end]*area[end])
#print("tau_conf = "); println(tau_conf)

tau_conf = sum(Vprime.*nash.*Tash)*(r[2]-r[1])
tau_conf = tau_conf / (hfluxash[end]*area[end])
#print("tau_conf = "); println(tau_conf)

#####
# Ash sim

dir = "ashsim"
f0alpha_ashsim = getdata("f0alpha")
f0ash_ashsim = getdata("f0ash")

ir = 1
plt[:figure]()
plot(v/vti[ir],vec(f0alpha_ashsim[ir,:]),"+b")
plot(v/vti[ir],vec(f0ash_ashsim[ir,:]),"-k")
xlabel(L"$v / v_{ti}$")
ylabel(L"Distribution function $F_{0}$ (A.U.)")
legend((L"$F_{0\alpha}$ from T3CORE",L"$F_\mathrm{ash}$"))
savefig("fit_lowen.png",bbox_inches="tight")
cla()
clf()
close()


########
# Dilution
#######

dir = "dilute_1e18_2"
heat_2 = vec(getdata("totheating"))
hflux_2 = vec(getdata("hflux"))
n_2 = vec(getdata("nalpha"))
T_2 = vec(getdata("Talpha"))

dir = "dilute_1e18_3"
heat_3 = vec(getdata("totheating"))
hflux_3 = vec(getdata("hflux"))
n_3 = vec(getdata("nalpha"))
T_3 = vec(getdata("Talpha"))


plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,heat_1e18*1.e-3,"-k")
plot(r/a,heat_2*1.e-3,"-.b")
plot(r/a,heat_3*1.e-3,"--m")
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(n_1e18.*T_1e18)*a/p0,"-k")
plot(r/a,getgradient(n_2.*T_2)*a/p0,"-.b")
plot(r/a,getgradient(n_3.*T_3)*a/p0,"--m")
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux_1e18*1.e-3,"-k")
plot(r/a,hflux_2*1.e-3,"-.b")
plot(r/a,hflux_3*1.e-3,"--m")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")
legend((L"$n_{\alpha,\mathrm{edge}} = n_e /100 $","Diluted source","Diluted source + turbulence"),fontsize=16,loc="lower right")

savefig("dilution_1e18.png",bbox_inches="tight")
cla()
clf()
close()


dir = "nedge_5e18"
heat_5e18 = vec(getdata("totheating"))
hflux_5e18 = vec(getdata("hflux"))
n_5e18 = vec(getdata("nalpha"))
T_5e18 = vec(getdata("Talpha"))

dir = "dilute_5e18_2"
heat_2 = vec(getdata("totheating"))
hflux_2 = vec(getdata("hflux"))
n_2 = vec(getdata("nalpha"))
T_2 = vec(getdata("Talpha"))

dir = "dilute_5e18_3"
heat_3 = vec(getdata("totheating"))
hflux_3 = vec(getdata("hflux"))
n_3 = vec(getdata("nalpha"))
T_3 = vec(getdata("Talpha"))


plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,heat_5e18*1.e-3,"-k")
plot(r/a,heat_2*1.e-3,"-.b")
plot(r/a,heat_3*1.e-3,"--m")
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(n_5e18.*T_5e18)*a/p0,"-k")
plot(r/a,getgradient(n_2.*T_2)*a/p0,"-.b")
plot(r/a,getgradient(n_3.*T_3)*a/p0,"--m")
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux_5e18*1.e-3,"-k")
plot(r/a,hflux_2*1.e-3,"-.b")
plot(r/a,hflux_3*1.e-3,"--m")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")
legend((L"$n_{\alpha,\mathrm{edge}} = n_e /20 $","Diluted source","Diluted source + turbulence"),fontsize=16,loc="lower right")

savefig("dilution_5e18.png",bbox_inches="tight")
cla()
clf()
close()

dir = "nedge_5e17"
heat_5e17 = vec(getdata("totheating"))
hflux_5e17 = vec(getdata("hflux"))
n_5e17 = vec(getdata("nalpha"))
T_5e17 = vec(getdata("Talpha"))

dir = "dilute_5e17_2"
heat_2 = vec(getdata("totheating"))
hflux_2 = vec(getdata("hflux"))
n_2 = vec(getdata("nalpha"))
T_2 = vec(getdata("Talpha"))

dir = "dilute_5e17_3"
heat_3 = vec(getdata("totheating"))
hflux_3 = vec(getdata("hflux"))
n_3 = vec(getdata("nalpha"))
T_3 = vec(getdata("Talpha"))


plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,heat_5e17*1.e-3,"-k")
plot(r/a,heat_2*1.e-3,"-.b")
plot(r/a,heat_3*1.e-3,"--m")
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(n_5e17.*T_5e17)*a/p0,"-k")
plot(r/a,getgradient(n_2.*T_2)*a/p0,"-.b")
plot(r/a,getgradient(n_3.*T_3)*a/p0,"--m")
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux_5e17*1.e-3,"-k")
plot(r/a,hflux_2*1.e-3,"-.b")
plot(r/a,hflux_3*1.e-3,"--m")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")
legend((L"$n_{\alpha,\mathrm{edge}} = n_e /200 $","Diluted source","Diluted source + turbulence"),fontsize=16,loc="lower right")

savefig("dilution_5e17.png",bbox_inches="tight")
cla()
clf()
close()

dir = "nedge_2e18"
heat_2e18 = vec(getdata("totheating"))
hflux_2e18 = vec(getdata("hflux"))
n_2e18 = vec(getdata("nalpha"))
T_2e18 = vec(getdata("Talpha"))

dir = "dilute_2e18_2"
heat_2 = vec(getdata("totheating"))
hflux_2 = vec(getdata("hflux"))
n_2 = vec(getdata("nalpha"))
T_2 = vec(getdata("Talpha"))

dir = "dilute_2e18_3"
heat_3 = vec(getdata("totheating"))
hflux_3 = vec(getdata("hflux"))
n_3 = vec(getdata("nalpha"))
T_3 = vec(getdata("Talpha"))


plt[:figure](figsize=(24,6))
plt[:tight_layout]
plt[:subplot](1,3,1)
plt[:tight_layout]
plot(r/a,heat_2e18*1.e-3,"-k")
plot(r/a,heat_2*1.e-3,"-.b")
plot(r/a,heat_3*1.e-3,"--m")
title("Heating",fontsize=24)
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,2)
plt[:tight_layout]
plot(r/a,getgradient(n_2e18.*T_2e18)*a/p0,"-k")
plot(r/a,getgradient(n_2.*T_2)*a/p0,"-.b")
plot(r/a,getgradient(n_3.*T_3)*a/p0,"--m")
title(L"$\nabla p_\alpha$",fontsize=24)
ylabel(L"$\alpha$ pressure gradient  $\, \times a / p_{e0}$",fontsize=16)
xlabel(L"$\psi / \psi_a $")

plt[:subplot](1,3,3)
plt[:tight_layout]
plot(r/a,hflux_2e18*1.e-3,"-k")
plot(r/a,hflux_2*1.e-3,"-.b")
plot(r/a,hflux_3*1.e-3,"--m")
title("Heat flux",fontsize=24)
ylabel(L"$\alpha$ heat flux $\left( \mathrm{kW}/\mathrm{m}^2 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")
legend((L"$n_{\alpha,\mathrm{edge}} = n_e /50 $","Diluted source","Diluted source + turbulence"),fontsize=16,loc="lower right")

savefig("dilution_2e18.png",bbox_inches="tight")
cla()
clf()
close()

dir = "splinek1"
ionheat_k1 = vec(getdata("ionheating"))

plt[:figure]()
plot(r/a,ionheating*1.e-3,"-k")
plot(r/a,ionheat_k1*1.e-3,"--b")
ylabel(L"Ion heating rate $\mathrm{kW} / \mathrm{m}^3$")
xlabel(L"$\psi / \psi_a $")
legend(("Cubic splines","Linear splines"))
savefig("splines.png",bbox_inches="tight")
cla()
clf()
close()


# Find "effective particle diffusion coefficients" for ash
Dr = Array(Float64,Nrad)
chi = Array(Float64,Nrad)
for ir in 1:Nrad
Dr[ir] = dot(d3v,vec(Drr[ir,:].*f0ash[ir,:]))/nash[ir]
chi[ir] = dot(d3v,vec(Drr[ir,:].*f0ash[ir,:]).*((energy/Tash[ir])-1.5).*energy)/(nash[ir]*Tash[ir])
end

plt[:figure]()
plot(r/a,Dr./chiGB)
plot(r/a,chi./chiGB)
ylabel(L"$D_{eff} /  \chi_{GB}$")
xlabel(L"$\psi / \psi_a $")
legend((L"$D_{ash,eff}$",L"$\chi_{ash,eff}$"))
savefig("Deff.png")
cla()
clf()
close()

plt[:tight_layout]
streamplot(r/a,v/valpha,rflux_eject',vflux_eject'*a/valpha)
xlabel(L"$\psi / \psi_a$",fontsize=36)
ylabel(L"$v / v_\alpha$",fontsize=36)
xlim(0.5,0.8)
ylim(0.0,1.05)
savefig("ejectstream.png",bbox_inches="tight")
cla()
clf()
close()

plt[:figure]()
plt[:tight_layout]
streamplot(r/a,v/valpha,rflux_weak',vflux_weak'*a/valpha)
#streamplot(r/a,v/valpha,rflux_weak',vflux_weak',color=sqrt(rflux_weak.^2+(vflux_weak).^2)')
#title(L"$\chi_{i,\mathrm{nominal}} / 5$"*"\n",fontsize=24)
#xlabel(L"$\psi / \psi_a$",fontsize=36)
#ylabel(L"$v / v_\alpha$",fontsize=36)
xlim(0.5,0.8)
ylim(0.0,1.05)
savefig("weakstream.png",bbox_inches="tight")
cla()
clf()
close()

plt[:tight_layout]
plot(r/a,elheating*1.e-3,"-k")
plot(r/a,elheat_strong*1.e-3,"--r")
plot(r/a,elheat_weak*1.e-3,":b")
ylabel(L"$\alpha$ heating $\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=16)
xlabel(L"$\psi / \psi_a $")
legend(("Nominal amplitude",L"$\chi_i \times 5$",L"$\chi_i / 5$"),fontsize=16)
savefig("elheatingstrength.png",bbox_inches="tight")
cla()
clf()
close()

PyCall.PyDict(matplotlib["rcParams"])["figure.subplot.hspace"] = "0.0"
plt[:figure](figsize=(6,6))
ax1 = plt[:subplot](2,1,1)
setp( ax1[:get_xticklabels](), visible=false)
plot(r/a,(sdiheating+sdeheating)*1.e-3,"--g")
plot(r/a,heat_weak*1.e-3,":b")
plot(r/a,totheating*1.e-3,"-k")
plot(r/a,heat_strong*1.e-3,"--r")
ylim(1.0,250.0)
xlim(0.5,0.8)
text(0.52,220.0,"a)",fontsize=18)
legend(("No turbulence",L"$\chi_i \times$ 0.2",L"Nominal $\chi_i$",L"$\chi_i \times$ 5"),loc="upper right",fontsize=14)
ylabel(L"$\alpha$ collisional heating"*"\n"*L"$\left( \mathrm{kW}/\mathrm{m}^3 \right)$",fontsize=14)

plt[:subplot](2,1,2)
plot(r/a,-getgradient(nsd.*Tsd)*a/p0,"--g")
plot(r/a,-getgradient(n_weak.*T_weak)*a/p0,":b")
plot(r/a,-getgradient(nalpha.*Talpha)*a/p0,"-k")
plot(r/a,-getgradient(n_strong.*T_strong)*a/p0,"--r")
#ylim(0.11,0.79)
ylim(0.11,1.0)
xlim(0.5,0.8)
text(0.52,0.85,"b)",fontsize=18)
ylabel("\n"*L"$ -\left(\partial p_\alpha / \partial r \right) \, \times a / p_{e0}$",fontsize=18)
xlabel(L"$\psi / \psi_a $",fontsize=18)

savefig("radialplots.png",bbox_inches="tight")
cla()
clf()
close()



