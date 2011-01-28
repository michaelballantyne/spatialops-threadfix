get_filename_component( dir ${CMAKE_CURRENT_LIST_FILE} PATH )

set( Particle_SRC
  ${dir}/ParticleFieldTypes.cpp
  )

if( ENABLE_TESTS )
  add_library( spatialops-particle
    STATIC
    ${Particle_SRC}
    )
  add_subdirectory( test )
endif( ENABLE_TESTS )

install( FILES
  ParticleFieldTypes.h
  ParticleOperators.h
  DESTINATION include/spatialops/particles
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  )

set( SO_SRC ${SO_SRC} ${Particle_SRC} PARENT_SCOPE )