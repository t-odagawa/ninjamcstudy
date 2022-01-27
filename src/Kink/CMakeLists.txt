message (STATUS "Kink study...")

add_executable(KinkSandMuon
	KinkSandMuon.cpp
	)

target_link_libraries(KinkSandMuon
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

# install the execute in the bin folder
install(TARGETS KinkSandMuon DESTINATION "${CMAKE_INSTALL_BINDIR}")
