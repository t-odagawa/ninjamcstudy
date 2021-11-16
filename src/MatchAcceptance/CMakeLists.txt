message (STATUS "Matching Acceptance...")

add_executable(MuonExpectedDistribution
	MuonExpectedDistribution.cpp
	)

target_link_libraries(MuonExpectedDistribution
	${ROOT_LIBRARIES}
	${Geant4_LIBRARIES}
	${Boost_LIBRARIES}
	${B2MC_LIBRARY}
	${NTBM_LIBRARY}
)

# install the execute in the bin folder
install(TARGETS MuonExpectedDistribution DESTINATION "${CMAKE_INSTALL_BINDIR}")