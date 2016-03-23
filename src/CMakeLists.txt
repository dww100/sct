cmake_minimum_required (VERSION 2.6)

#SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
#SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX})

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
enable_language (Fortran)

add_executable(brktos brktos.f90 sjp_util.f90)
add_executable(brksph2pdb brksph2pdb.f90 sjp_util.f90)
add_executable(rfacXN rfacXN2.f90 sjp_util.f90)
add_executable(aps aps.f90)
add_executable(frcc frcc.f90)
add_executable(hypro hypro.f90)
add_executable(sct sct.f90)
add_executable(sluv sluv.f90)
add_executable(sctpl sctpl.f axtrans.f decide.f getfile.f grinit.f linegt.f opnfil.f plots.f ptypget.f rfctxt.f smear.f toalph.f widang.f consgt.f errmsg.f getfval.f guinlab.f linfit.f plotone.f prtscr.f qcnvt.f sctpl.f symbgt.f unaxtr.f datype.f fitq.f getmodel.f init.f lotus.f plotq.f pstyle.f rfacauto.f setaxes.f titles.f wghtdav.f)

SET (EXES brktos brksph2pdb rfacXN aps frcc hypro sct sluv sctpl)
SET (PERMS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(TARGETS ${EXES} PERMISSIONS ${PERMS} DESTINATION . COMPONENT binaries)
