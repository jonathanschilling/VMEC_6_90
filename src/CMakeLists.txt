
add_subdirectory(coilslib)
add_subdirectory(ga_opt)
add_subdirectory(miscel)
add_subdirectory(lapack)
add_subdirectory(vmec2000)
add_subdirectory(vmec_auxmod)

set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
