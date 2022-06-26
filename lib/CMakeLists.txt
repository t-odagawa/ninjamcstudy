add_library(NinjaMcStudy SHARED
			 HistogramStyle.hpp
			 HistogramStyle.cpp
  			 MyFunction.hpp
  			 MyFunction.cpp
			 DrawConst.hpp
)

target_include_directories(NinjaMcStudy PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB NINJAMC_LIB_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
set_target_properties(NinjaMcStudy PROPERTIES PUBLIC_HEADER "${NINJAMC_LIB_INCLUDES}")

install(TARGETS NinjaMcStudy
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCDIR}")