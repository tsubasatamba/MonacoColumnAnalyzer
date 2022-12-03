include_directories(
  include
  ${ANLG4_INC_DIR}
  ${G4_INC_DIR}
  ${CFITSIO_INC_DIR}
  ${ROOT_INC_DIR}
  $ENV{HOME}/include/BinaryAnalyzer
  )

link_directories(
  /usr/local/lib
  ${ANLG4_LIB_DIR}
  )

set(MCA_CORE_CLASSES
  src/MaterialProperties.cc
  src/EventCollection.cc
  src/EmissionProfile.cc
  )

add_library(MCACore SHARED
  ${MCA_CORE_CLASSES}
  )

target_link_libraries(MCACore
  ${ANLG4_LIB}
  ${G4_LIB}
  ${CFITSIO_LIB}
  ${ROOT_LIB}
  $ENV{HOME}/lib/libBACore.dylib
  )

install(TARGETS MCACore LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/MonacoColumnAnalyzer
  FILES_MATCHING PATTERN "*.hh")
