
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/reconstruct.f
	${CMAKE_CURRENT_SOURCE_DIR}/splines.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmodules.f
	${CMAKE_CURRENT_SOURCE_DIR}/vvacuum.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
