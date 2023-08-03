project(shad-cpp0)

cmake_minimum_required(VERSION 3.8)

set(CMAKE_CXX_STANDARD             20)
set(CMAKE_MODULE_PATH              "${CMAKE_SOURCE_DIR}/tools/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_EXPORT_COMPILE_COMMANDS  ON)

find_package(Catch REQUIRED)
find_package(PNG)
find_package(JPEG)

find_package(Poco QUIET COMPONENTS Foundation Net JSON)
if (NOT Poco_FOUND)
  message(STATUS "Seems like POCO is not installed on your machine.")
  message(STATUS "Tasks that depend on POCO are disabled.")
endif()

find_package(PythonLibs 3 QUIET)
if (NOT PythonLibs_FOUND)
  message(STATUS "python3 libraries not found. Python tasks are disabled.")
endif()

include(tools/cmake/TestSolution.cmake)
include(tools/cmake/BuildFlags.cmake)
include_directories(tools/util)

function(add_task name)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${name})
    add_subdirectory(${name})
  endif()
endfunction()

add_task(multiplication)

add_task(bin-pow)
add_task(range)
add_task(unique)
add_task(split)
add_task(polish-notation)
add_task(word-count)
add_task(permutations)
add_task(reverse-map)
add_task(long-sum)
add_task(diff-pairs)

add_task(quadratic)
add_task(pt-triangle)
add_task(factorization)
add_task(rotate)
add_task(sort-students)
add_task(war)
add_task(filter-even)
add_task(swap-sort)
add_task(entrance)

add_task(decoder)

add_task(ring-buffer)
add_task(matrix)
add_task(static-map)
add_task(string-view)
add_task(dungeon)
add_task(stack)

add_task(bad-hash)

add_task(find-last)
add_task(local-max)
add_task(partition)
add_task(fold)
add_task(functors)
add_task(itertools)

if (${PNG_FOUND} AND ${JPEG_FOUND})
  add_task(raytracer-debug)
  add_task(raytracer-geom)
  add_task(raytracer-reader)
  add_task(raytracer)
  add_task(raytracer-b2)
elseif()
  message(STATUS "raytracer disabled. PNG:${PNG_FOUND} JPEG:${JPEG_FOUND}")
endif()

add_task(bad-rand)

add_task(vector)
add_task(lru-cache)
add_task(cow-vector)

add_task(dedup)
add_task(smart-ptr)
add_task(deque)
add_task(rule-of-5)
add_task(gcd)
add_task(reallol)

add_task(list)
add_task(calculator)
add_task(reader)

add_task(scheme)
add_task(scheme-tokenizer)
add_task(scheme-parser)

add_task(safe-transform)
add_task(intrusive-list)
add_task(any)

add_task(tryhard)
add_task(pimpl)

add_task(editor)
add_task(small-test-framework)

if (Poco_FOUND)
  add_task(hello-http)
  add_task(bot)
endif()
