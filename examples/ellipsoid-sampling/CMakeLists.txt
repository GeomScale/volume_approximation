# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis
# Copyright (c) 2021- Vaibhav Thakkar

# Contributed and/or modified by Vaibhav Thakkar

# Licensed under GNU LGPL.3, see LICENCE file


add_executable (ellipsoid2d-sampling ellipsoid2d-sampling.cpp)
TARGET_LINK_LIBRARIES(ellipsoid2d-sampling ${LP_SOLVE})