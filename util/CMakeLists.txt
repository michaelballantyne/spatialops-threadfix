
add_library( timelogger STATIC TimeLogger.cpp )
target_link_libraries( timelogger ${Boost_LIBRARIES} )
install( TARGETS timelogger DESTINATION lib/spatialops )
install( FILES TimeLogger.h DESTINATION include/spatialops/test )

