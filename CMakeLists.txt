cmake_minimum_required(VERSION 3.11)
project(sim_indirect_recip_solitary)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
#set(CMAKE_VERBOSE_MAKEFILE 1)


# See http://google.github.io/googletest/quickstart-cmake.html
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
enable_testing()

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_MakeAvailable(json)

# Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)

add_executable(test_Norm test_Norm.cpp Norm.hpp)
target_link_libraries(test_Norm GTest::gtest_main)
add_test(NAME test_Norm COMMAND test_Norm)

add_executable(test_PrivRepGame test_PrivRepGame.cpp PrivRepGame.hpp Norm.hpp)
target_link_libraries(test_PrivRepGame GTest::gtest_main)
add_test(NAME test_PrivRepGame COMMAND test_PrivRepGame)

add_executable(test_SolitaryObsGame test_SolitaryObsGame.cpp SolitaryObsGame.hpp Norm.hpp)
target_link_libraries(test_SolitaryObsGame GTest::gtest_main)
add_test(NAME test_SolitaryObsGame COMMAND test_SolitaryObsGame)

add_executable(inspect_Norm inspect_Norm.cpp Norm.hpp)

add_executable(inspect_PrivRepGame inspect_PrivRepGame.cpp PrivRepGame.hpp Norm.hpp)
target_link_libraries(inspect_PrivRepGame nlohmann_json::nlohmann_json)

add_executable(inspect_HscoreDynamics inspect_HscoreDynamics.cpp SolitaryObsGame.hpp Norm.hpp)
target_link_libraries(inspect_HscoreDynamics nlohmann_json::nlohmann_json)

add_executable(inspect_SolitaryObsGame inspect_SolitaryObsGame.cpp SolitaryObsGame.hpp Norm.hpp)
target_link_libraries(inspect_SolitaryObsGame nlohmann_json::nlohmann_json)

add_executable(main_CorrelatedOpinions main_CorrelatedOpinions.cpp PrivRepGame.hpp Norm.hpp)
target_link_libraries(main_CorrelatedOpinions nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

add_executable(main_PayoffComparison main_PayoffComparison.cpp PrivRepGame.hpp Norm.hpp)
target_link_libraries(main_PayoffComparison nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

add_executable(inspect_PrivRepGameTernary inspect_PrivRepGameTernary.cpp PrivRepGameTernary.hpp NormTernary.hpp)
target_link_libraries(inspect_PrivRepGameTernary nlohmann_json::nlohmann_json)

add_executable(inspect_PrivRepGame3rd inspect_PrivRepGame3rd.cpp PrivRepGame3rd.hpp Norm3rd.hpp)
target_link_libraries(inspect_PrivRepGame3rd nlohmann_json::nlohmann_json)
