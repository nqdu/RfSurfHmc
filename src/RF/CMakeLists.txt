
file(GLOB SOURCE "*.cpp" "*.f90")
pybind11_add_module(librf ${SOURCE})

target_link_libraries(librf PRIVATE ${FFTW_LIBRARIES} )
install(TARGETS librf
  COMPONENT python
  LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/model/lib
)
# add_library(libcps330 ${SOURCE})