# CMakeList.txt: проект CMake для LabaNM3; включите исходный код и определения,
# укажите здесь логику для конкретного проекта.
#
cmake_minimum_required (VERSION 3.8)

set (TARGET_NAME LabaNM3)

# Добавьте источник в исполняемый файл этого проекта.
add_executable (${TARGET_NAME} "LabaNM3.cpp")

# set binary filenames
set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME_DEBUG ${TARGET_NAME}-debug )
set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME_RELWITHDEBINFO ${TARGET_NAME}-relwidthdebinfo )
set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME_RELEASE ${TARGET_NAME}-release )
set_target_properties( ${TARGET_NAME} PROPERTIES OUTPUT_NAME_MINSIZEREL ${TARGET_NAME}-minsizerel )

# set warning levels
if (MSVC)
	target_compile_options(${TARGET_NAME} PUBLIC /W4 /permissive)
else()
	target_compile_options(${TARGET_NAME} PUBLIC -Wall -Wextra -pedantic)
endif()