ttk_add_vtk_module()

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(ttkEireneWriter PRIVATE TTK_ENABLE_EIGEN)
  target_include_directories(ttkEireneWriter SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})
endif()
