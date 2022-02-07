
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/coilslib.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmodules_coils.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
