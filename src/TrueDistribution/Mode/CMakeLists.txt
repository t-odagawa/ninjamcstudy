message (STATUS "True distribution by mode...")

add_executable(MuonTrueDistributionMode
	MuonTrueDistributionMode.cpp
	)

target_link_libraries(MuonTrueDistributionMode
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

add_executable(ProtonTrueDistributionMode
	ProtonTrueDistributionMode.cpp
)

target_link_libraries(ProtonTrueDistributionMode
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

add_executable(PionTrueDistributionMode
	PionTrueDistributionMode.cpp
)

target_link_libraries(PionTrueDistributionMode
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

# install the execute in the bin folder
install(TARGETS MuonTrueDistributionMode DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS ProtonTrueDistributionMode DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(TARGETS PionTrueDistributionMode DESTINATION "${CMAKE_INSTALL_BINDIR}")
