ttk_add_base_library(implicitTriangulation
  SOURCES
    ImplicitTriangulation.cpp
    ImplicitPreconditions.cpp
  HEADERS
    ImplicitTriangulation.h
    ImplicitPreconditions.h
  DEPENDS
    abstractTriangulation
    periodicImplicitTriangulation
    )

    # target_compile_options(implicitTriangulation PRIVATE -fsanitize=address -O0)
    # target_link_options(implicitTriangulation PRIVATE -fsanitize=address)
