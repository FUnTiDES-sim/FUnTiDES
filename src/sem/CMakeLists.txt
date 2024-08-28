################################
# Setup for SEM_EXECUTABLE
################################

if(NOT DEFINED DEVICE)
  set (SEM_EXECUTABLE "sem_${SOLVER}.exe" CACHE PATH "" FORCE)
else()
  set (SEM_EXECUTABLE "sem_${SOLVER}_${DEVICE}.exe" CACHE PATH "" FORCE)
endif()

message ( "\n---------------------------------------------------------------------------------" )
message ( "-- The Current Program Model Used is: ${SOLVER}")
message ( "-- The Current executable to be built is: ${SEM_EXECUTABLE}")
message ( "-- CMAKE_SOURCE_DIR : ${CMAKE_SOURCE_DIR}")
message ( "-- PROJECT_BINARY_DIR :  ${PROJECT_BINARY_DIR}")
message ( "---------------------------------------------------------------------------------\n" )

#-------------------------------------------------------
# SEM proxy application setup
#-------------------------------------------------------

set (sem_sources
     finiteElements/SEMQkGL.cpp
     mesh/SEMmesh.cpp
     solver/SEMsolver.cpp
     main/SEMproxy.cpp
     main/main.cpp)

blt_add_executable( NAME ${SEM_EXECUTABLE}
                    SOURCES ${sem_sources}
                    DEPENDS_ON ${extraComponentsLinkList} ${dependencyList} )

target_include_directories( ${SEM_EXECUTABLE} PUBLIC
        ${CMAKE_BINARY_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/../utils
        ${CMAKE_CURRENT_SOURCE_DIR}/finiteElements
        ${CMAKE_CURRENT_SOURCE_DIR}/mesh
        ${CMAKE_CURRENT_SOURCE_DIR}/solver
        ${CMAKE_CURRENT_SOURCE_DIR}/main)

install(TARGETS ${SEM_EXECUTABLE}
        RUNTIME DESTINATION bin)

