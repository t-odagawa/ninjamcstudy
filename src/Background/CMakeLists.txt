message (STATUS "MC BG distribution...")

# Distribution
add_executable(ProtonModuleBG
	ProtonModuleBG.cpp
	)

target_link_libraries(ProtonModuleBG
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS ProtonModuleBG DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(WallBG
	WallBG.cpp
	)

target_link_libraries(WallBG
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS WallBG DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(BabyMindBG
	BabyMindBG.cpp
	)

target_link_libraries(BabyMindBG
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS BabyMindBG DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(DownWagasciBG
	DownWagasciBG.cpp
	)

target_link_libraries(DownWagasciBG
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS DownWagasciBG DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(OtherEccBG
	OtherEccBG.cpp
	)

target_link_libraries(OtherEccBG
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS OtherEccBG DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(IronIntBG
	IronIntBG.cpp
	)

target_link_libraries(IronIntBG
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS IronIntBG DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(AnuIntBG
	AnuIntBG.cpp
	)

target_link_libraries(AnuIntBG
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS AnuIntBG DESTINATION "${CMAKE_INSTALL_BINDIR}")

message(STATUS "Momentum Consistency...")

add_executable(MomentumConsistency
	MomentumConsistency.cpp
	)

target_link_libraries(MomentumConsistency
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS MomentumConsistency DESTINATION "${CMAKE_INSTALL_BINDIR}")
