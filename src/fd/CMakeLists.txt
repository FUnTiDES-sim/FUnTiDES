################################
# Setup for FD_EXECUTABLE
################################

if(NOT DEFINED DEVICE)
  set (FD_EXECUTABLE "fd_${SOLVER}.exe" CACHE PATH "" FORCE)
else()
  set (FD_EXECUTABLE "fd_${SOLVER}_${DEVICE}.exe" CACHE PATH "" FORCE)
endif()

message ( "\n---------------------------------------------------------------------------------" )
message ( "-- The Current Program Model Used is: ${SOLVER}")
message ( "-- The Current executable to be built is: ${FD_EXECUTABLE}")
message ( "-- CMAKE_CURRENT_SOURCE_DIR :${CMAKE_CURRENT_SOURCE_DIR}")
message ( "-- CMAKE_BINARY_DIR :${CMAKE_BINARY_DIR}")
message ( "---------------------------------------------------------------------------------\n" )

#-------------------------------------------------------
# fd proxy application setup
#-------------------------------------------------------

#set (fd_sources FDTD${SOLVER}.cpp) 
set (fd_sources FDTDmain.cpp) 

blt_add_executable( NAME ${FD_EXECUTABLE}
                    SOURCES ${fd_sources}
                    DEPENDS_ON ${extraComponentsLinkList} ${dependencyList} )

target_include_directories( ${FD_EXECUTABLE} PUBLIC
        ${CMAKE_BINARY_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/utilsFDTD 
        ${CMAKE_CURRENT_SOURCE_DIR}/../utils )

install(TARGETS ${FD_EXECUTABLE}
        RUNTIME DESTINATION bin)
