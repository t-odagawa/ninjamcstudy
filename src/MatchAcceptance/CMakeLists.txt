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

install(TARGETS MuonExpectedDistribution DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(TangentDistributionPaper
	TangentDistributionPaper.cpp
	)

target_link_libraries(TangentDistributionPaper
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS TangentDistributionPaper DESTINATION "${CMAKE_INSTALL_BINDIR}")