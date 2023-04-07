cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(SealPIR VERSION 2.1 LANGUAGES CXX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(main
	sha256.cpp
)

# target_link_libraries()

add_library(sealpir STATIC
  pir.cpp
  pir_client.cpp
  pir_server.cpp
)

find_package(SEAL 3.2.0 EXACT REQUIRED)
#find_package(SEAL 4.0 EXACT REQUIRED)        // edit by ZL 20230214, since the version of SEAL is 4.0 now

target_link_libraries(main sealpir SEAL::seal)
target_link_libraries(main gmp pbc ssl crypto)
