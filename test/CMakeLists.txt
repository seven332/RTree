include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE
  GIT_SUBMODULES "")
FetchContent_MakeAvailable(googletest)

add_executable(rtree_test rtree_test.cpp)
target_link_libraries(rtree_test rtree gtest_main gmock)
target_include_directories(rtree_test PRIVATE . ../src)

include(GoogleTest)
gtest_discover_tests(rtree_test)
