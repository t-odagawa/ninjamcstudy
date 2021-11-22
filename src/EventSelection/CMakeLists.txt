message (STATUS "Event selection study...")

add_executable(MinimumDistance
	MinimumDistance.cpp
	MinimumDistance.hpp
	)

target_link_libraries(MinimumDistance
	${ROOT_LIBRARIES}
	${Geant4_LIBRARIES}
	${Boost_LIBRARIES}
	${B2MC_LIBRARY}
	${NTBM_LIBRARY}
)

# install the execute in the bin folder
install(TARGETS MinimumDistance DESTINATION "${CMAKE_INSTALL_BINDIR}")