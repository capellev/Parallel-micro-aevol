# Copyright (C) 2018  INRIA
#
# Treerecs is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Treerecs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Build GTEST but not GMOCK, install neither
option(BUILD_GTEST "Builds the googletest subproject" ON)
option(BUILD_GMOCK "Builds the googlemock subproject" OFF)
option(INSTALL_GTEST "Disable installation of googletest" OFF)
option(INSTALL_GMOCK "Disable installation of googlemock" OFF)
add_subdirectory(googletest)

# Libraries to link unit tests against
set(test_libs
  gtest_main
    mini_aevol)

# List unit tests
set(TESTS
    IndividualTest)

# Create a runner for each unit test generate TEST_RUNNERS
foreach (TEST IN LISTS TESTS)
  set(TEST_RUNNER  run_${TEST})
  set(TEST_RUNNERS ${TEST_RUNNERS} ${TEST_RUNNER})
  add_executable(${TEST} ${TEST}.cpp)
  set_target_properties(${TEST}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
    )
  target_link_libraries(${TEST} PUBLIC ${test_libs})
  add_custom_target(${TEST_RUNNER}
    COMMAND ${CMAKE_BINARY_DIR}/tests/${TEST}
    DEPENDS ${TEST})
  set_target_properties(${TEST_RUNNER}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
    )
endforeach(TEST)

# Add target check
add_custom_target(check DEPENDS utest)

# Add target utest
add_custom_target(utest DEPENDS ${TEST_RUNNERS})

cuda_add_executable(threefrytest threefry.cu ../Threefry.cpp)
target_link_libraries(threefrytest ${ZLIB_LIBRARIES} )
