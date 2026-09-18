add_executable(test_LastKrylovSlots test_LastKrylovSlots.cpp)
target_link_libraries(test_LastKrylovSlots PRIVATE Catch2::Catch2WithMain psimaglite)
target_include_directories(test_LastKrylovSlots PRIVATE ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(test_LastKrylovSlots)

add_executable(testVectorWithOffsets testVectorWithOffsets.cpp)
target_link_libraries(testVectorWithOffsets PRIVATE Catch2::Catch2WithMain dmrgpp_utils psimaglite)
target_include_directories(testVectorWithOffsets PRIVATE ${CMAKE_SOURCE_DIR}/dmrg ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(testVectorWithOffsets)

add_executable(test_PermutationParity test_PermutationParity.cpp)
target_link_libraries(test_PermutationParity PRIVATE Catch2::Catch2WithMain psimaglite)
target_include_directories(test_PermutationParity PRIVATE ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(test_PermutationParity)

# BatchedGemm test
add_executable(test_BatchedGemm test_BatchedGemm.cpp)
target_link_libraries(test_BatchedGemm PRIVATE Catch2::Catch2WithMain dmrgpp_utils psimaglite)
if(PLUGIN_SC)
  target_link_libraries(test_BatchedGemm PRIVATE gpuplugin)
endif()
target_include_directories(test_BatchedGemm PRIVATE ${CMAKE_SOURCE_DIR}/dmrg/Engine
                                                    ${CMAKE_SOURCE_DIR}/dmrg/Engine/MatrixVectorKron)
catch_discover_tests(test_BatchedGemm)

add_executable(test_MatrixVectorStored test_MatrixVectorStored.cpp)
target_link_libraries(test_MatrixVectorStored PRIVATE Catch2::Catch2WithMain dmrgpp_utils psimaglite)
target_include_directories(test_MatrixVectorStored PRIVATE ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(test_MatrixVectorStored)

add_executable(bench_BatchedGemm bench_BatchedGemm.cpp)
target_link_libraries(bench_BatchedGemm PRIVATE dmrgpp_utils psimaglite)
target_include_directories(bench_BatchedGemm PRIVATE ${CMAKE_SOURCE_DIR}/dmrg/Engine
                                                     ${CMAKE_SOURCE_DIR}/dmrg/Engine/MatrixVectorKron)
