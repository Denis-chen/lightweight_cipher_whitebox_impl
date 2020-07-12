
#1.cmake verson，指定cmake版本
cmake_minimum_required(VERSION 3.2)

#2.project name，指定项目的名称，一般和项目的文件夹名称对应
PROJECT(WHITE_PAPER)

#3.head file path，头文件目录
INCLUDE_DIRECTORIES(include)


#4.source directory，源文件目录
AUX_SOURCE_DIRECTORY(library LIB_SRCS)

# AUX_SOURCE_DIRECTORY(library/sms4 LIB_SRCS)
SET(LIB_SRCS ${LIB_SRCS} library/sms4/sms4_common.c library/sms4/sms4_setkey.c library/sms4/sms4_enc.c)

# test
SET(TEST_EXEC ${TEST_EXEC} test/utils/testutil.c)
SET(TEST_EXEC ${TEST_EXEC} test/sm4_control_experiment.cpp test/NTLTest.cpp test/klein_vector.cpp test/LBlock_test.cpp test/present_test.cpp)

SET(LIB_SRCS ${LIB_SRCS} library/affine/NTLUtils.cpp )
SET(LIB_SRCS ${LIB_SRCS} library/affine/MatrixUtils.cpp )


#klein
AUX_SOURCE_DIRECTORY(library/klein LIB_SRCS)

#LBlock
AUX_SOURCE_DIRECTORY(library/LBlock LIB_SRCS)

#present
SET(LIB_SRCS
    ${LIB_SRCS}
    library/present/present80.cpp
    library/present/present80wb.cpp    
)
# AUX_SOURCE_DIRECTORY(library/present LIB_SRCS)


AUX_SOURCE_DIRECTORY(library/internal LIB_SRCS)
AUX_SOURCE_DIRECTORY(library/sms4wb LIB_SRCS)

#AUX_SOURCE_DIRECTORY(test TEST_EXEC)

AUX_SOURCE_DIRECTORY(. TEST_EXEC)


# set name for build
SET(BUILD_NAME ${PROJECT_NAME})
SET(BUILD_TEST_NAME ${PROJECT_NAME}_TEST)
SET(TEST_EXEC_NAME ${PROJECT_NAME}_MAIN)



SET(CXXFLAGS $ENV{CXXFLAGS})

#set(CMAKE_CXX_FLAGS "-std=c11 -stdlib=libc++")

#add_compile_options(-std=c11)

SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -std=c++11 -O0")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")


# SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -Os -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer/SDKs/iPhoneOS10.3.sdk -F/Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer/SDKs/iPhoneOS10.3.sdk/System/Library/Frameworks -arch arm64 -arch armv7 -miphoneos-version-min=8.0 -flto")

ADD_DEFINITIONS(-DDUMMY_ROUND)
ADD_DEFINITIONS(-Ddebug)

#GMP
SET(GMP_DIRS "/usr/local/lib/libgmp.dylib")

# NTL
if (NTL_NO_DEFAULT_PATH)
    set(NTL_MODIF NO_DEFAULT_PATH)
    message("Not using default paths for NTL search, ${NTL_MODIF}")
endif (NTL_NO_DEFAULT_PATH)
message("HOME: $ENV{HOME}")
find_path(NTL_INCLUDE_PATH NAMES NTL/mat_GF2.h PATHS $ENV{HOME}/ntl $ENV{HOME}/ntl/include /usr/local /opt/local /usr ${NTL_MODIF})
message("NTL_INCLUDE_PATH: ${NTL_INCLUDE_PATH}")

find_library(NTL_LIB ntl PATHS ${CMAKE_CURRENT_SOURCE_DIR} $ENV{HOME}/ntl $ENV{HOME}/ntl/lib /usr/local /opt/local /usr ${NTL_MODIF})
if(NOT NTL_LIB)
    message(FATAL_ERROR "ntl library not found.  Rerun cmake with -DCMAKE_PREFIX_PATH=\"<path to lib1>;<path to lib2>\"")
endif()
message("NTL_LIB: ${NTL_LIB}")
message("source files: ${LIB_SRCS}")

include_directories(${NTL_INCLUDE_PATH})


# Linking
# target_link_libraries(main ${NTL_LIB} ${GMP_DIRS} ${Boost_LIBRARIES})
# target_link_libraries(testing ${NTL_LIB} ${GMP_DIRS} ${Boost_LIBRARIES})



#6.add executable file，添加要编译的可执行文件
ADD_EXECUTABLE(${TEST_EXEC_NAME} ${LIB_SRCS} ${TEST_EXEC}  )
TARGET_LINK_LIBRARIES(${TEST_EXEC_NAME} ${NTL_LIB} ${GMP_DIRS} ${Boost_LIBRARIES})

