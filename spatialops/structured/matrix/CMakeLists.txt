if( ENABLE_TESTS )
  add_subdirectory( test )
endif( ENABLE_TESTS )


# only install if this is enabled.
if( ENABLE_MATRIX )

  install(
    FILES
    FVRestrictOp.h
    FVStaggeredBCOp.h
    FVStaggeredDivergence.h
    FVStaggeredGradient.h
    FVStaggeredInterpolant.h
    FVStaggeredIndexHelper.h
    FVStaggeredScratch.h
    FVTopHatFilter.h
    DESTINATION include/spatialops/structured/matrix/
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
    )

  install( FILES
    FVStaggered.h
    FVStaggeredOperatorTypes.h
    DESTINATION include/spatialops/structured/
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
    )

endif( ENABLE_MATRIX )
