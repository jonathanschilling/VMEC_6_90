
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/vmodules.f
	${CMAKE_CURRENT_SOURCE_DIR}/gade_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/kind_spec.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
