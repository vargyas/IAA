cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(IAA)
#Module
set(MODULE IAA)
add_definitions(-D_MODULE_="${MODULE}")


list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS}
                    ${CMAKE_SOURCE_DIRS}
  )
add_definitions(${ROOT_CXX_FLAGS})

# Sources - alphabetical order
set(SRCS
  AliJHistManagerROOT6.cxx
  )

# Headers from sources
set(HDRS
  AliJHistManagerROOT6.h
    JFiete.h
    loadFilip.h
    loadFilipIAAfinal.h
    mcorr.h
    mfit.h
    miaa.h
    mplot.h
    mtools.h
)


# Generate the dictionary
ROOT_GENERATE_DICTIONARY(G__${MODULE} ${HDRS} LINKDEF IAALinkDef.h)

# Add a shared library
add_library(${MODULE} SHARED   G__${MODULE}.cxx)
#add_library(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Linking the library
target_link_libraries(${MODULE} ${ROOT_LIBRARIES})

add_executable(process processDraw.cxx)
target_link_libraries(process ${MODULE}) 



# Public include folders that will be propagated to the dependecies
#target_include_directories(${MODULE} PUBLIC ${incdirs})

# System dependent: Modify the way the library is build
#if(${CMAKE_SYSTEM} MATCHES Darwin)
#	set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
#endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
#	install(TARGETS ${MODULE}
        #			ARCHIVE DESTINATION lib
            #			LIBRARY DESTINATION lib)
        #install(FILES ${HDRS} DESTINATION include)
