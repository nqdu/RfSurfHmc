
file(GLOB SOURCE "*.cpp" "*.f90" "*.f")
pybind11_add_module(libsurf ${SOURCE})

target_link_libraries(libsurf PRIVATE ${FFTW_LIBRARIES} )
install(TARGETS libsurf
  COMPONENT python
  LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/model/lib
)
# add_library(libcps330 ${SOURCE})