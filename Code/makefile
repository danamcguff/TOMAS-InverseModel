# Makefile for TOMAS box model
# questions or comments, email Jeff Pierce, jrpierce@andrew.cmu.edu
###########################################################################

# options
comp = gfortran
# to run quickly, use -O2 flag (optimize level 2) or -O if using g77 with no -g flag (for debugger)
# to debug, remove -O2 and add -g flag
o_flags = -c -g -fbounds-check -fcheck=all# flags for the compiling of *.o files
exe_flags = -g -fbounds-check # flags for the compiling of *.exe files

ROOTDIR := /usr/lib
#PARDISO := /Users2/danamcg/FORTRANfunctions
LinAlg := -L$(ROOTDIR) -lblas -llapack
#LinAlg2 := -L$(PARDISO)-lpardiso500-GNU481-X86-64


# groups
box_obj  = array.o	box.o	boxinit.o	getlen.o	\
	   boxchem.o	setoh.o		emisnaclarke_box.o      \
	   logemiss.o   soacond.o      select_inventory.o        \
	   feedback.o affine_dynamics.o remaining_dynamics.o 
cond_obj = tmcond.o 	dmdt_int.o	nucleation.o	\
	   aqoxid.o	vehk_nucl.o	napa_nucl.o	getGrowthTime.o \
	   getCoagLoss.o		cond_nuc.o	ezcond.o	\
	   getCondSink.o		getH2SO4prod.o	getNucConc.o	\
	   getNucRate.o	getH2SO4conc.o	getNucBinFit.o	nucCond.o	\
	   NH3_GISStoTOMAS.o
coag_obj = multicoag.o
size_obj = initbounds.o	loginit.o	report.o	speciesmap.o	\
	   storenm.o	aerodiag.o	ezwatereqm.o	eznh3eqm.o	\
	   mnfix.o	scalemom.o	waterso4.o	waternacl.o
prop_obj = gasdiff.o	aerodens.o


# executables
box.exe : $(box_obj) $(cond_obj) $(coag_obj) $(size_obj) $(prop_obj) BB192SM9.COM BT263box.COM sizecode.COM makefile
	$(comp)  $(exe_flags) -o box.exe $(box_obj) $(cond_obj) $(coag_obj) $(size_obj) $(prop_obj) $(emis_obj) -llapack 


# object files
array.o:array.f BB192SM9.COM BT263box.COM sizecode.COM makefile
	$(comp) $(o_flags) array.f
aerodiag.o:aerodiag.f BB192SM9.COM BT263box.COM sizecode.COM makefile
	$(comp) $(o_flags) aerodiag.f
aqoxid.o : aqoxid.f BB192SM9.COM BT263box.COM sizecode.COM makefile
	$(comp) $(o_flags) aqoxid.f
box.o : box.f BB192SM9.COM BT263box.COM sizecode.COM  makefile
	$(comp) $(o_flags) box.f
boxinit.o : boxinit.f BB192SM9.COM BT263box.COM makefile
	$(comp) $(o_flags) boxinit.f
boxchem.o : boxchem.f makefile
	$(comp) $(o_flags) boxchem.f
setoh.o : setoh.f makefile
	$(comp) $(o_flags) setoh.f
tmcond.o : tmcond.f sizecode.COM makefile
	$(comp) $(o_flags) tmcond.f
getlen.o : getlen.f makefile
	$(comp) $(o_flags) getlen.f
nucleation.o : nucleation.f BB192SM9.COM BT263box.COM sizecode.COM makefile
	$(comp) $(o_flags) nucleation.f
emisnaclarke_box.o : emisnaclarke_box.f BB192SM9.COM BT263box.COM sizecode.COM makefile
	$(comp) $(o_flags) emisnaclarke_box.f
multicoag.o : multicoag.f sizecode.COM makefile
	$(comp) $(o_flags) multicoag.f
dmdt_int.o : dmdt_int.f makefile
	$(comp) $(o_flags) dmdt_int.f
initbounds.o : initbounds.f sizecode.COM makefile
	$(comp) $(o_flags) initbounds.f
logemiss.o : logemiss.f BB192SM9.COM BT263box.COM sizecode.COM makefile
	$(comp) $(o_flags) logemiss.f
loginit.o : loginit.f BB192SM9.COM BT263box.COM sizecode.COM makefile
	$(comp) $(o_flags) loginit.f
report.o : report.f BB192SM9.COM BT263box.COM makefile
	$(comp) $(o_flags) report.f
speciesmap.o : speciesmap.f BB192SM9.COM BT263box.COM makefile
	$(comp) $(o_flags) speciesmap.f
storenm.o : storenm.f sizecode.COM makefile
	$(comp) $(o_flags) storenm.f
ezwatereqm.o : ezwatereqm.f sizecode.COM makefile
	$(comp) $(o_flags) ezwatereqm.f
eznh3eqm.o : eznh3eqm.f sizecode.COM makefile
	$(comp) $(o_flags) eznh3eqm.f
gasdiff.o  : gasdiff.f makefile
	$(comp) $(o_flags) gasdiff.f
aerodens.o : aerodens.f makefile
	$(comp) $(o_flags) aerodens.f
mnfix.o : mnfix.f sizecode.COM makefile
	$(comp) $(o_flags) mnfix.f
scalemom.o : scalemom.f BB192SM9.COM BT263box.COM makefile
	$(comp) $(o_flags) scalemom.f
waterso4.o : waterso4.f makefile
	$(comp) $(o_flags) waterso4.f
waternacl.o : waternacl.f makefile
	$(comp) $(o_flags) waternacl.f
napa_nucl.o : napa_nucl.f makefile
	$(comp) $(o_flags) napa_nucl.f
vehk_nucl.o : vehk_nucl.f makefile
	$(comp) $(o_flags) vehk_nucl.f
getGrowthTime.o : getGrowthTime.f makefile
	$(comp) $(o_flags) getGrowthTime.f
getCoagLoss.o : getCoagLoss.f makefile
	$(comp) $(o_flags) getCoagLoss.f
cond_nuc.o : cond_nuc.f makefile
	$(comp) $(o_flags) cond_nuc.f
ezcond.o : ezcond.f makefile
	$(comp) $(o_flags) ezcond.f
getCondSink.o : getCondSink.f makefile
	$(comp) $(o_flags) getCondSink.f
getH2SO4prod.o : getH2SO4prod.f makefile
	$(comp) $(o_flags) getH2SO4prod.f
getNucConc.o : getNucConc.f makefile
	$(comp) $(o_flags) getNucConc.f
getNucRate.o : getNucRate.f makefile
	$(comp) $(o_flags) getNucRate.f
getH2SO4conc.o : getH2SO4conc.f makefile
	$(comp) $(o_flags) getH2SO4conc.f
nucCond.o : nucCond.f makefile
	$(comp) $(o_flags) nucCond.f
getNucBinFit.o : getNucBinFit.f makefile
	$(comp) $(o_flags) getNucBinFit.f
NH3_GISStoTOMAS.o : NH3_GISStoTOMAS.f makefile
	$(comp) $(o_flags) NH3_GISStoTOMAS.f
feedback.o : feedback.f makefile
	$(comp) $(o_flags) feedback.f
select_inventory.o : select_inventory.f makefile
	$(comp) $(o_flags) select_inventory.f
affine_dynamics.o : affine_dynamics.f makefile
	$(comp) $(o_flags) affine_dynamics.f
remaining_dynamics_inventory.o : remaining_dynamics_inventory.f makefile
	$(comp) $(o_flags) remaining_dynamics_inventory.f
remaining_dynamics.o : remaining_dynamics.f makefile
	$(comp) $(o_flags) remaining_dynamics.f
soacond.o : soacond.f makefile
	$(comp) $(o_flags) soacond.f
clean :
	rm so2_* coagsink_so2_* condsink_so2_* parameter_so2_* NucRate_so2_* measurements_so2_* soacond_so2_* pbio_data* pbio_debug* SizeDist_so2*
