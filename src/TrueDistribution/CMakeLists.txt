message (STATUS "True distribution...")

add_executable(MuonTrueDistribution
	MuonTrueDistribution.cpp
	)

target_link_libraries(MuonTrueDistribution
	${ROOT_LIBRARIES}
	${Geant4_LIBRARIES}
	${Boost_LIBRARIES}
	${B2MC_LIBRARY}
)

add_executable(ProtonTrueDistribution
	ProtonTrueDistribution.cpp
	)

target_link_libraries(ProtonTrueDistribution
	${ROOT_LIBRARIES}
	${Geant4_LIBRARIES}
	${Boost_LIBRARIES}
	${B2MC_LIBRARY}
)

add_executable(PionTrueDistribution
	PionTrueDistribution.cpp
	)

target_link_libraries(PionTrueDistribution
	${ROOT_LIBRARIES}
	${Geant4_LIBRARIES}
	${Boost_LIBRARIES}
	${B2MC_LIBRARY}
)

# install the execute in the bin folder
install(TARGETS MuonTrueDistribution DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS ProtonTrueDistribution DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS PionTrueDistribution DESTINATION "${CMAKE_INSTALL_BINDIR}")
