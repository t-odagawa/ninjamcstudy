message (STATUS "MC distribution...")

# Distribution
add_executable(MomchDraw
	MomchDraw.cpp
	AnalyzeMultiplicity.hpp
	AnalyzeMultiplicity.cpp
	AnalyzeVertex.hpp
	AnalyzeVertex.cpp
	AnalyzeMomentum.hpp
	AnalyzeMomentum.cpp
	AnalyzeAngle.hpp
	AnalyzeAngle.cpp
	AnalyzePid.hpp
	AnalyzePid.cpp
	)

target_link_libraries(MomchDraw
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS MomchDraw DESTINATION "${CMAKE_INSTALL_BINDIR}")

message (STATUS "Angle correction factor...")

add_executable(AngleCorrection
	AngleCorrection.cpp)

target_link_libraries(AngleCorrection
	${NINJA_MCSTUDY_LIBRARIES}
	NinjaMcStudy
)

install(TARGETS AngleCorrection DESTINATION "${CMAKE_INSTALL_BINDIR}")