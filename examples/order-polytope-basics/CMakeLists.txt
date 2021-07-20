# VolEsti (volume computation and sampling library)

add_executable (order_polytope order_polytope.cpp)
TARGET_LINK_LIBRARIES(order_polytope ${LP_SOLVE})
