project(UnscentedKF)

cmake_minimum_required (VERSION 3.5)


#FIND_PACKAGE(PythonLibs 2.7 REQUIRED)
#MESSAGE("${PYTHON_LIBRARY}")
#MESSAGE("${PYTHON_INCLUDE_DIR}")

add_definitions(-std=c++0x)
include_directories(/usr/include/python2.7/)

set(sources
   src/ukf.cpp
   src/main.cpp )

add_executable(UnscentedKF ${sources} )
#target_link_libraries(UnscentedKF /usr/local/Cellar/python/2.7.13/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib)