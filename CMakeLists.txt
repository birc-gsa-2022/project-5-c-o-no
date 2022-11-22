cmake_minimum_required (VERSION 3.21)
project (PROJECT)
enable_testing()

file(GLOB readmapFiles
        "src/func/*.c"
        "src/func/parsers/*.c"
        )

add_executable(readmap src/readmap.c)
add_test(example example)

add_executable(testing1 test/tests.c ${readmapFiles})
add_executable(testing2 test/testSearch.c ${readmapFiles})
add_executable(testRotate test/testRotation.c ${readmapFiles})
add_executable(time test/time.c ${readmapFiles})
add_test(tester testing1 "./")
add_test(tester2 testing2 "./")
add_test(testRotate testRotate "./")
