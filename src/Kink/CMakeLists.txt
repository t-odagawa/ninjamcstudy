message (STATUS "Kink study...")

# Kink sandmuon
add_executable(KinkSandMuon
	KinkSandMuon.cpp
	)

target_link_libraries(KinkSandMuon
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS KinkSandMuon DESTINATION "${CMAKE_INSTALL_BINDIR}")

# Kink interaction
add_executable(KinkInteraction
	KinkInteraction.cpp
	)

target_link_libraries(KinkInteraction
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS KinkSandMuon DESTINATION "${CMAKE_INSTALL_BINDIR}")

# Kink distribution
add_executable(MakeKinkDistribution
	MakeKinkDistribution.cpp
	)

target_link_libraries(MakeKinkDistribution
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS MakeKinkDistribution DESTINATION "${CMAKE_INSTALL_BINDIR}")

# Momch merge
add_executable(MomchMerge
	MomchMerge.cpp
	MomchMerge.hpp
	)

target_link_libraries(MomchMerge
	${NINJA_MCSTUDY_LIBRARIES}
)

install(TARGETS MomchMerge DESTINATION "${CMAKE_INSTALL_BINDIR}")
