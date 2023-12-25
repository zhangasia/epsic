set(SOURCES 
    test_Dirac.cpp     
    test_inner_product.cpp   
    test_Matrix.cpp      
    test_myfinite.cpp     
    test_rotation.cpp 
    test_Basis.cpp         
    test_eigen.cpp      
    test_Jacobi.cpp          
    test_Minkowski.cpp   
    test_Pauli.cpp        
    test_Stokes.cpp 
    test_Conventions.cpp   
    test_Estimate.cpp   
    test_Jones.cpp           
    test_Mueller.cpp     
    test_Quaternion.cpp   
    test_Vector.cpp 
)
foreach(SOURCE ${SOURCES})
    get_filename_component(NAME ${SOURCE} NAME_WE)
    # message(${SOURCE})
    add_executable(${NAME} ${SOURCE})
    target_link_libraries(${NAME} PUBLIC epsiclib)
endforeach()
