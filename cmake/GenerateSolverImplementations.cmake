# ============================================================================
# GenerateSolverImplementations.cmake
# ============================================================================
# Function to generate solver implementation files from templates
# ============================================================================

function(generate_solver_implementations)
    set(options "")
    set(oneValueArgs OUTPUT_DIR OUTPUT_VAR SOLVER_NAME FILE_PREFIX)
    set(multiValueArgs ORDERS MESH_TYPES MODEL_TYPES PHYSIC_TYPES)
    cmake_parse_arguments(GEN "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    # Validate required arguments
    if(NOT DEFINED GEN_SOLVER_NAME)
        message(FATAL_ERROR "SOLVER_NAME argument is required")
    endif()
    if(NOT DEFINED GEN_ORDERS)
        message(FATAL_ERROR "ORDERS argument is required")
    endif()
    if(NOT DEFINED GEN_MESH_TYPES)
        message(FATAL_ERROR "MESH_TYPES argument is required")
    endif()
    if(NOT DEFINED GEN_MODEL_TYPES)
        message(FATAL_ERROR "MODEL_TYPES argument is required")
    endif()
    if(NOT DEFINED GEN_PHYSIC_TYPES)
        message(FATAL_ERROR "PHYSIC_TYPES argument is required")
    endif()
    if(NOT DEFINED GEN_OUTPUT_DIR)
        set(GEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
    endif()
    if(NOT DEFINED GEN_OUTPUT_VAR)
        message(FATAL_ERROR "OUTPUT_VAR argument is required")
    endif()
    if(NOT DEFINED GEN_FILE_PREFIX)
        set(GEN_FILE_PREFIX "SEMQ")
    endif()

    # Initialize output list
    set(GENERATED_SOURCES)
    
    # Generate solver implementation files
    foreach(ORDER ${GEN_ORDERS})
        foreach(MESH_TYPE ${GEN_MESH_TYPES})
            foreach(PHYSIC_TYPE ${GEN_PHYSIC_TYPES})
                set(TEMPLATE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/templates/${GEN_SOLVER_NAME}_solver_template_${PHYSIC_TYPE}.cpp.in")
                
                foreach(MODEL_TYPE ${GEN_MODEL_TYPES})
                    set(GEN_FILENAME "${GEN_FILE_PREFIX}${ORDER}_${MESH_TYPE}_${GEN_SOLVER_NAME}_solver_${PHYSIC_TYPE}_model${MODEL_TYPE}.cpp")
                    set(GEN_FILEPATH "${GEN_OUTPUT_DIR}/${GEN_FILENAME}")
                    
                    # Determine model header and class based on mesh type
                    if(MESH_TYPE STREQUAL "struct")
                        set(MODEL_HEADER "model_struct.h")
                        set(MODEL_CLASS "ModelStruct")
                    else()
                        set(MODEL_HEADER "model_unstruct.h")
                        set(MODEL_CLASS "ModelUnstruct")
                    endif()
                    
                    # Determine if this is "onelements" (false) or "onnodes" (true) for last template param
                    if(MODEL_TYPE STREQUAL "onelements")
                        set(MODEL_ON_NODES "false")
                    else()
                        set(MODEL_ON_NODES "true")
                    endif()
                    
                    # Determine if ModelStruct needs ORDER template parameter
                    if(MESH_TYPE STREQUAL "struct")
                        set(MODEL_TEMPLATE_PARAMS "float, int, ORDER")
                    else()
                        set(MODEL_TEMPLATE_PARAMS "float, int")
                    endif()
                    
                    # Configure the file from template
                    configure_file(
                        "${TEMPLATE_FILE}"
                        "${GEN_FILEPATH}"
                        @ONLY
                    )
                    
                    # Add to list of sources
                    list(APPEND GENERATED_SOURCES "${GEN_FILEPATH}")
                    
                    message(STATUS "Configured: ${GEN_FILENAME}")
                endforeach()
            endforeach()
        endforeach()
    endforeach()
    
    # Return the list of generated sources to the parent scope
    set(${GEN_OUTPUT_VAR} ${GENERATED_SOURCES} PARENT_SCOPE)
endfunction()

# ============================================================================
# Function to generate p-adaptive solver implementation files from templates
# ============================================================================
# Unlike generate_solver_implementations(), the p-adaptive solver couples two
# sub-domains of different polynomial order (ORDER_MIN, ORDER_MAX), so one
# file is generated per ordered pair ORDER_MIN < ORDER_MAX drawn from ORDERS.
# ============================================================================

function(generate_padaptive_solver_implementations)
    set(options "")
    set(oneValueArgs OUTPUT_DIR OUTPUT_VAR SOLVER_NAME FILE_PREFIX)
    set(multiValueArgs ORDERS MESH_TYPES MODEL_TYPES PHYSIC_TYPES)
    cmake_parse_arguments(GEN "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # Validate required arguments
    if(NOT DEFINED GEN_SOLVER_NAME)
        message(FATAL_ERROR "SOLVER_NAME argument is required")
    endif()
    if(NOT DEFINED GEN_ORDERS)
        message(FATAL_ERROR "ORDERS argument is required")
    endif()
    if(NOT DEFINED GEN_MESH_TYPES)
        message(FATAL_ERROR "MESH_TYPES argument is required")
    endif()
    if(NOT DEFINED GEN_MODEL_TYPES)
        message(FATAL_ERROR "MODEL_TYPES argument is required")
    endif()
    if(NOT DEFINED GEN_PHYSIC_TYPES)
        message(FATAL_ERROR "PHYSIC_TYPES argument is required")
    endif()
    if(NOT DEFINED GEN_OUTPUT_DIR)
        set(GEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
    endif()
    if(NOT DEFINED GEN_OUTPUT_VAR)
        message(FATAL_ERROR "OUTPUT_VAR argument is required")
    endif()
    if(NOT DEFINED GEN_FILE_PREFIX)
        set(GEN_FILE_PREFIX "PADAPTIVEQ")
    endif()

    # Initialize output list
    set(GENERATED_SOURCES)

    # Generate solver implementation files for every ORDER_MIN < ORDER_MAX pair
    foreach(ORDER_MIN ${GEN_ORDERS})
        foreach(ORDER_MAX ${GEN_ORDERS})
            if(ORDER_MIN LESS ORDER_MAX)
                foreach(MESH_TYPE ${GEN_MESH_TYPES})
                    foreach(PHYSIC_TYPE ${GEN_PHYSIC_TYPES})
                        set(TEMPLATE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/templates/${GEN_SOLVER_NAME}_solver_template_${PHYSIC_TYPE}.cpp.in")

                        foreach(MODEL_TYPE ${GEN_MODEL_TYPES})
                            set(GEN_FILENAME "${GEN_FILE_PREFIX}${ORDER_MIN}_${ORDER_MAX}_${MESH_TYPE}_${GEN_SOLVER_NAME}_solver_${PHYSIC_TYPE}_model${MODEL_TYPE}.cpp")
                            set(GEN_FILEPATH "${GEN_OUTPUT_DIR}/${GEN_FILENAME}")

                            # Determine model header and class based on mesh type
                            if(MESH_TYPE STREQUAL "struct")
                                set(MODEL_HEADER "model_struct.h")
                                set(MODEL_CLASS "ModelStruct")
                            else()
                                set(MODEL_HEADER "model_unstruct.h")
                                set(MODEL_CLASS "ModelUnstruct")
                            endif()

                            # Determine if this is "onelements" (false) or "onnodes" (true) for last template param
                            if(MODEL_TYPE STREQUAL "onelements")
                                set(MODEL_ON_NODES "false")
                            else()
                                set(MODEL_ON_NODES "true")
                            endif()

                            # ModelStruct is built on the pMax mesh (highest order); ModelUnstruct is order-agnostic
                            if(MESH_TYPE STREQUAL "struct")
                                set(MODEL_TEMPLATE_PARAMS "float, int, ORDER_MAX")
                            else()
                                set(MODEL_TEMPLATE_PARAMS "float, int")
                            endif()

                            # Configure the file from template
                            configure_file(
                                "${TEMPLATE_FILE}"
                                "${GEN_FILEPATH}"
                                @ONLY
                            )

                            # Add to list of sources
                            list(APPEND GENERATED_SOURCES "${GEN_FILEPATH}")

                            message(STATUS "Configured: ${GEN_FILENAME}")
                        endforeach()
                    endforeach()
                endforeach()
            endif()
        endforeach()
    endforeach()

    # Return the list of generated sources to the parent scope
    set(${GEN_OUTPUT_VAR} ${GENERATED_SOURCES} PARENT_SCOPE)
endfunction()
