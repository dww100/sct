cmake_minimum_required(VERSION 2.6)

find_program(PYTHON "python")

if (PYTHON)

    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/sct/__init__.py")
    set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    add_custom_command(OUTPUT ${OUTPUT}
    COMMAND ${PYTHON} ${SETUP_PY} build
    COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
    DEPENDS ${DEPS})
    add_custom_target(target ALL DEPENDS ${OUTPUT})
    if (CMAKE_INSTALL_PREFIX)
        install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install --user --install-scripts=${CMAKE_INSTALL_PREFIX})")
    else()
        install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")
    endif()

endif()
