# include the standard dependencies
INCLUDE_DIRECTORIES(SYSTEM ${ALLPIX_DEPS_INCLUDE_DIRS})

# define the library
ADD_LIBRARY(AllpixObjects SHARED
    ChargeDeposit.cpp
)

# link the standard dependencies
TARGET_LINK_LIBRARIES(AllpixObjects ${ALLPIX_LIBRARIES} ${ALLPIX_DEPS_LIBRARIES})

# specify install for the messages
INSTALL(TARGETS AllpixObjects
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib)