using input:input_for_analytic_test, Nrad, Nv, a, deltat, D0, rgrid_in
using grids:init_vgrid,init_rgrid, rgrid, v
using geometry:init_geometry
using matrix: zero_collop,build_matrix, f0, gindex, solve_steadystate, global_matrix, advance_timestep, init_collop
using sourcemod: init_source_for_analytic_test
using diffcoeff: init_diffcoeff_for_analytic_test
using species: init_species_for_analytic_test
using boundary: calculate_boundary_for_analytic_test, F0edge
using constants: valpha
using PyPlot

#plt.style[:use]("myfig")
#plt[:style][:use]("myfig")

println("Performing analytic test...")
input_for_analytic_test()
init_rgrid()
init_vgrid()
init_geometry()
#init_species_for_tests()
zero_collop()
init_source_for_analytic_test()
init_diffcoeff_for_analytic_test()
build_matrix()
calculate_boundary_for_analytic_test()

if deltat > 0.0

  nt = 30
  for it in 1:nt
    println("Timestep ",it)
    advance_timestep(it)
  end

else

  solve_steadystate()

end

f0anal = zeros(Nrad,Nv)
f0pred = zeros(Nrad,Nv)
relerr = zeros(Nrad,Nv)
abserr = zeros(Nrad,Nv)
for ir in 1:Nrad
  for iv in 1:Nv
    idx = gindex(ir,iv)
    f0anal[ir,iv] = (2.0-rgrid[ir]^2/a^2)*exp(-v[iv]^2/valpha^2)
    f0pred[ir,iv] = f0[idx] #+ F0edge[iv]
    relerr[ir,iv] = abs(f0anal[ir,iv]-f0pred[ir,iv])/abs(f0anal[ir,iv])
    abserr[ir,iv] = abs(f0anal[ir,iv]-f0pred[ir,iv])
  end
end

# r-norm:
# For each velocity, sum error over all radii, then take maximum
errVSv = zeros(Nv)
for iv in 1:Nv
   errVSv[iv] = mean(abserr[:,iv])
end
rnorm = maximum(errVSv)

# v-norm:
# For each radius, sum error over all velocities, then take maximum
errVSr = zeros(Nrad)
for ir in 1:Nrad
   errVSr[ir] = mean(vec(abserr[ir,1:end]))
end
vnorm = maximum(errVSr)

absnorm = norm(abserr)
relnorm = norm(relerr)

#println("Abs. norm = ",absnorm)
#println("Rel. norm = ",relnorm)
println("Max abs error = ",maximum(abserr))
println("Max rel error = ",maximum(relerr[:,1:end-1]))

f0compfile = open("f0comp.dat","w+")
#writedlm(f0file, [v/valpha f0anal[1,:]' f0anal[end/2,:]'])
for iv in 1:Nv
  writedlm(f0compfile,[v[iv]/valpha,f0pred[1,iv],f0pred[end/2,iv]]')
end
close(f0compfile)

resfile = open("resolutions.dat","a")
writedlm(resfile,[Nrad,Nv,rnorm,vnorm,norm(abserr),norm(relerr),maximum(abserr),maximum(relerr),norm(abserr[2:end-1,2:end-1]),norm(relerr[2:end-1,2:end-1]),maximum(abserr[2:end-1,2:end-1]),maximum(relerr[2:end-1,2:end-1])]')
close(resfile)

pcolormesh(rgrid[2:end-1]/a,v[2:end-1]/valpha,abserr[2:end-1,2:end-1]')
xlabel(L"$r / a$")
ylabel(L"$v / v_\mathrm{ref}$")
colorbar()
savefig("internalabserrs.png")
clf()
cla()

plot(v/valpha,vec(f0pred[1,:]),"+r")
plot(v/valpha,vec(f0anal[1,:]),"-k",lw=2)
xlabel(L"$v / v_\alpha$")
ylabel(L"$f_0$")
legend(("Prediction","Analytic"),fontsize=12)
savefig("r0.png")
clf()
cla()

plt[:tight_layout]
plt[:subplots_adjust](wspace=0.4,bottom=0.2)
subplot(1,2,1)
plot(v[1:end/2]/valpha,vec(f0pred[1,1:end/2]),":+k")
plot(v[1:end/2]/valpha,vec(f0pred[3,1:end/2]),":+r")
plot(v[1:end/2]/valpha,vec(f0pred[5,1:end/2]),":+g")
plot(v[1:end/2]/valpha,vec(f0pred[8,1:end/2]),":+b")
plot(v[1:end/2]/valpha,vec(f0anal[1,1:end/2]),"-k",linewidth=1)
plot(v[1:end/2]/valpha,vec(f0anal[3,1:end/2]),"-r",linewidth=1)
plot(v[1:end/2]/valpha,vec(f0anal[5,1:end/2]),"-g",linewidth=1)
plot(v[1:end/2]/valpha,vec(f0anal[8,1:end/2]),"-b",linewidth=1)
xlabel(L"$v / v_\mathrm{ref}$")
ylabel(L"$F_0$")
#legend(("Prediction","Analytic"),fontsize=12)
legend(("ir=1","ir=3","ir=5","ir=8"),fontsize=12)
#savefig("rcomp.png")

subplot(1,2,2)
plot(rgrid/a,f0pred[:,1],":+k")
plot(rgrid/a,f0pred[:,indmin(abs(valpha-v))],":+g")
plot(rgrid/a,f0pred[:,indmin(abs(2*valpha-v))],":+r")
plot(rgrid/a,f0anal[:,1],"-k",linewidth=1)
plot(rgrid/a,f0anal[:,indmin(abs(valpha-v))],"-g",linewidth=1)
plot(rgrid/a,f0anal[:,indmin(abs(2*valpha-v))],"-r",linewidth=1)
legend(("iv = 1","iv = "string(indmin(abs(valpha-v))) ,"iv = "string(indmin(abs(2*valpha-v)))),fontsize=12)
xlabel(L"$r / a$")
ylabel(L"$F_0$")
#legend(("Prediction","Analytic"),fontsize=12)
savefig("analytictest_rv.png")
clf()
cla()

plt[:tight_layout]
semilogy(v[1:end]/valpha,vec(f0pred[1,1:end]),":+k")
semilogy(v[1:end]/valpha,vec(f0pred[8,1:end]),":+b")
semilogy(v[1:end]/valpha,vec(f0anal[1,1:end]),"-k",linewidth=1)
semilogy(v[1:end]/valpha,vec(f0anal[8,1:end]),"-b",linewidth=1)
xlabel(L"$v / v_\mathrm{ref}$")
ylabel(L"$F_0$")
xlim(0.0,5.0)
legend(("ir=1","ir=8"))
savefig("vcomp_log.png")

pcolormesh(rgrid/a,v[1:Nv-1]/valpha,relerr[:,1:Nv-1]')
xlabel(L"$r / a$")
ylabel(L"$v / v_\mathrm{ref}$")
colorbar()
savefig("relerrs.png")
clf()
cla()

pcolormesh(rgrid/a,v/valpha,abserr')
xlabel(L"$r / a$")
ylabel(L"$v / v_\mathrm{ref}$")
xlim(0.3,1.0)
ylim(0.0,5.0)
colorbar()
savefig("abserrs.png")
clf()
cla()

