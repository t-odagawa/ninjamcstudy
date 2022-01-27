message (STATUS "Event selection study...")

add_executable(MinimumDistance
	MinimumDistance.cpp
	)

target_link_libraries(MinimumDistance
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

# install the execute in the bin folder
install(TARGETS MinimumDistance DESTINATION "${CMAKE_INSTALL_BINDIR}")