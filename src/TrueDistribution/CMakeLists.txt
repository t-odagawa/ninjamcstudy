message (STATUS "True distribution...")

add_executable(MuonTrueDistribution
	MuonTrueDistribution.cpp
	)

target_link_libraries(MuonTrueDistribution
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

add_executable(ProtonTrueDistribution
	ProtonTrueDistribution.cpp
	)

target_link_libraries(ProtonTrueDistribution
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

add_executable(PionTrueDistribution
	PionTrueDistribution.cpp
	)

target_link_libraries(PionTrueDistribution
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

# install the execute in the bin folder
install(TARGETS MuonTrueDistribution DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS ProtonTrueDistribution DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS PionTrueDistribution DESTINATION "${CMAKE_INSTALL_BINDIR}")
