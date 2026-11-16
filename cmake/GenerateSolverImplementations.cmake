# ============================================================================
# GenerateSolverImplementations.cmake
# ============================================================================
# Function to generate solver implementation files from templates
# ============================================================================

# Convenience wrapper for the most common instantiation pattern: one physics,
# all orders from 1 to a variable holding the maximum, both mesh types, both
# model locations, added to one library target. Each directory that used to
# carry a hand-written copy of this boilerplate now calls:
#
#   add_solver_instantiations(
#     PHYSIC      acoustic|elastic|acoustoelastic
#     MAX_ORDER   <variable holding the maximum polynomial order>
#     [SOLVER_NAME makutu|gemm]  integral back-end (default makutu)
#     [FILE_PREFIX ...]          file prefix override (default per generator)
#     [TARGET fun_solver]        target receiving the sources (default)
#     [ORDER_PAIRS]              p-adaptive: one file per ordered pair
#                                ORDER_MIN < ORDER_MAX instead of one per
#                                order (dispatches to the p-adaptive
#                                generator; requires the *_solver_template_*.cpp.in
#                                of the calling directory to use @ORDER_MIN@/
#                                @ORDER_MAX@)
#   )
function(add_solver_instantiations)
    set(options ORDER_PAIRS)
    set(oneValueArgs PHYSIC MAX_ORDER SOLVER_NAME FILE_PREFIX TARGET)
    cmake_parse_arguments(ASI "${options}" "${oneValueArgs}" "" ${ARGN})

    if(NOT DEFINED ASI_PHYSIC)
        message(FATAL_ERROR "add_solver_instantiations: PHYSIC is required")
    endif()
    if(NOT DEFINED ASI_MAX_ORDER)
        message(FATAL_ERROR "add_solver_instantiations: MAX_ORDER is required")
    endif()
    if(NOT DEFINED ASI_SOLVER_NAME)
        set(ASI_SOLVER_NAME "makutu")
    endif()
    if(NOT DEFINED ASI_TARGET)
        set(ASI_TARGET "fun_solver")
    endif()

    set(ORDERS "")
    foreach(order RANGE 1 ${${ASI_MAX_ORDER}})
        list(APPEND ORDERS ${order})
    endforeach()

    set(EXTRA_ARGS "")
    if(DEFINED ASI_FILE_PREFIX)
        list(APPEND EXTRA_ARGS FILE_PREFIX "${ASI_FILE_PREFIX}")
    endif()

    # Arguments common to both generators; the ORDER_PAIRS option only
    # selects which one iterates: generate_solver_implementations() emits one
    # file per order, generate_padaptive_solver_implementations() one per
    # ordered pair ORDER_MIN < ORDER_MAX.
    set(GEN_ARGS
        SOLVER_NAME "${ASI_SOLVER_NAME}"
        ${EXTRA_ARGS}
        ORDERS ${ORDERS}
        MESH_TYPES "struct" "unstruct"
        MODEL_TYPES "onelements" "onnodes"
        PHYSIC_TYPES "${ASI_PHYSIC}"
        OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}"
        OUTPUT_VAR ASI_GENERATED_SOURCES
    )
    if(ASI_ORDER_PAIRS)
        generate_padaptive_solver_implementations(${GEN_ARGS})
    else()
        generate_solver_implementations(${GEN_ARGS})
    endif()
    target_sources(${ASI_TARGET} PRIVATE ${ASI_GENERATED_SOURCES})
endfunction()

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
