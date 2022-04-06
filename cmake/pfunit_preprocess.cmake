#
# Given a path to a file return its name and it extension.
#
# For example given "/path/to/file.f90"
# Will set filename="file" and extension=".f90"
#
function(split_filepath filename extension path)

  # OBTAIN EXTENSION
  string(REGEX MATCH "\\.[a-zA-Z0-9]+$" ext ${path})

  # REMOVE EXTENSION
  string(REGEX REPLACE "${ext}$" "" path_temp ${path})

  # OBTAIN FILE NAME WITHOUT EXTENSION
  string(REGEX MATCH "[a-zA-Z0-9_]+$" name ${path_temp})

  set(${filename} ${name} PARENT_SCOPE)
  set(${extension} ${ext} PARENT_SCOPE)

endfunction(split_filepath)


#
# Perform preprocessing of unit tests
#
# Each file in a 'test_list' will be preprocessed by the pfUnit preprocessor and stored in
# 'folder_name'. In additions test suite registry "testSuites.inc" wil be created in 'folder_name'
#
# Arguments:
#   processed_list [out] : List of test source files after preprocessing
#   folder_name [in] : Name of the folder in which preprocessed files will be collected
#   test_list [in] : List of files to preprocessed
#
function(pfunit_preprocess processed_list folder_name test_list)

  # Create directory in binary_dir for temporary preprocessed test files
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${folder_name})

  # Create required test suites listing file
  file(WRITE ${PROJECT_BINARY_DIR}/${folder_name}/testSuites.inc "")

  # Preprocess collected test files from global property into variable
  # Copy files to the build folder
  set(tests)

  foreach(_testPath IN LISTS ${test_list})
    split_filepath(testName extension ${_testPath})

    # ADD RULE PREPROCESS ALL TEST FILES TO A SINGLE FOLDER IN BINARY_DIR
    add_custom_command(
       OUTPUT ${PROJECT_BINARY_DIR}/${folder_name}/${testName}${extension}
       COMMAND ${PYTHON_EXECUTABLE} ${PFUNIT_PREPROC}/pFUnitParser.py ${_testPath}
                                    ${PROJECT_BINARY_DIR}/${folder_name}/${testName}${extension}
       DEPENDS pfunit ${_testPath}
       COMMENT "Preprocessing test ${testName}"
       VERBATIM
       )
    # APPEND LIST OF ALL PREPROCESSED UNIT TEST FILES & ADD TEST SUITE TO testSuites.inc
    set(tests ${tests} ${PROJECT_BINARY_DIR}/${folder_name}/${testName}${extension})
    file(APPEND ${PROJECT_BINARY_DIR}/${folder_name}/testSuites.inc "ADD_TEST_SUITE(${testName}_suite)\n")
  endforeach()

  # Return variables to the parent scope
  set(${processed_list} ${tests} PARENT_SCOPE)

endfunction(pfunit_preprocess)
