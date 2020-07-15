cmake_minimum_required(VERSION 3.0.0)
project(FastGreeksWithAlgorithmicDifferentiation VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "-std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -llapack -lblas")

set(SOURCES "${CMAKE_SOURCE_DIR}/src")
set(INCLUDES "${CMAKE_SOURCE_DIR}/include")
set(DOC "${CMAKE_SOURCE_DIR}/doc")
set(BIN "${CMAKE_SOURCE_DIR}/bin")
set(DATA "${CMAKE_SOURCE_DIR}/data")
set(LIBRARIES "${CMAKE_SOURCE_DIR}/libraries")

set(EXECUTABLE_OUTPUT_PATH ${BIN})

include_directories(${INCLUDES} ${CMAKE_SYSTEM_INCLUDE_PATH})
include_directories(${LIBRARIES} ${CMAKE_SYSTEM_INCLUDE_PATH})

add_executable(demoTangentMode "${SOURCES}/demoTangentMode.cpp")
add_executable(demoAdjointMode "${SOURCES}/demoAdjointMode.cpp")

add_executable(basketOptionTime "${SOURCES}/basketOptionTime.cpp")
add_executable(basketOptionShow "${SOURCES}/basketOptionShow.cpp")

add_executable(bestOfAsianOptionTime "${SOURCES}/bestOfAsianOptionTime.cpp")
add_executable(bestOfAsianOptionShow "${SOURCES}/bestOfAsianOptionShow.cpp")

add_executable(correlationTime "${SOURCES}/correlationTime.cpp")
add_executable(correlationShow "${SOURCES}/correlationShow.cpp")


