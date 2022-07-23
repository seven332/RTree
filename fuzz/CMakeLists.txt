add_executable(rtree_fuzz rtree_fuzz.cpp)
target_compile_options(
  rtree_fuzz PRIVATE -fsanitize=address,undefined,fuzzer
                     -fprofile-instr-generate -fcoverage-mapping)
target_link_libraries(rtree_fuzz rtree)
target_link_options(rtree_fuzz PRIVATE -fsanitize=address,undefined,fuzzer
                    -fprofile-instr-generate -fcoverage-mapping)
