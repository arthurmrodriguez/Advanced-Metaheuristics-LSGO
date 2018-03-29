add_custom_target (sprng
		  COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/SRC > /dev/null 2> /dev/null
		  COMMENT "Building SPRNG library. Output has been redirected to /dev/null")
