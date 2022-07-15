message (STATUS "Efficiency details...")

# Distribution
add_executable(EventEfficiencyBreakdown
	EventEfficiencyBreakdown.cpp
	)

target_link_libraries(EventEfficiencyBreakdown
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS EventEfficiencyBreakdown DESTINATION "${CMAKE_INSTALL_BINDIR}")