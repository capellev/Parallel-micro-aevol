cmake_minimum_required(VERSION 3.8.2)
project(pdc_mini_aevol)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(SFMT-src-1.4)

add_definitions(-DSFMT_MEXP=607)

#find_package(CUDA REQUIRED)
#enable_language(CUDA)
#set(CUDA_NVCC_FLAGS "-std=c++11 -arch=compute_60 -lineinfo")

#cuda_compile(cuda_aevol Algorithms.cu ThreefryGPU.cu)

add_library(mini_aevol_gpu ExpManager.cpp ExpManager.h AeTime.cpp AeTime.h JumpingMT.cpp JumpingMT.h DnaMutator.cpp DnaMutator.h MutationEvent.cpp MutationEvent.cpp Organism.cpp Organism.h Promoter.h RNA.h Protein.h Stats.cpp Stats.h Gaussian.h Threefry.cpp Dna.cpp main.cpp)
target_include_directories(mini_aevol_gpu INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})


find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    #target_link_libraries(mini_aevol PUBLIC sfmt PUBLIC ${ZLIB_LIBRARIES})
    target_link_libraries(mini_aevol_gpu PUBLIC sfmt PUBLIC ${ZLIB_LIBRARIES})
endif( ZLIB_FOUND )

add_executable(pdc_mini_aevol main.cpp )
#cuda_add_executable(pdc_mini_aevol_gpu main.cpp ${cuda_aevol})
# nvToolsExt for enhanced profiling (ad-hoc chunks)
target_link_libraries(pdc_mini_aevol sfmt ${ZLIB_LIBRARIES} mini_aevol_gpu)

#target_link_libraries(pdc_mini_aevol_gpu PUBLIC mini_aevol_gpu)

#add_subdirectory(test EXCLUDE_FROM_ALL)
