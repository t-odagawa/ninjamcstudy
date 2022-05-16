message (STATUS "Kink study...")

# Kink sandmuon
add_executable(MomchDraw
	MomchDraw.cpp
	AnalyzeMultiplicity.hpp
	AnalyzeMultiplicity.cpp
	AnalyzeVertex.hpp
	AnalyzeVertex.cpp
	)

target_link_libraries(MomchDraw
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS MomchDraw DESTINATION "${CMAKE_INSTALL_BINDIR}")

