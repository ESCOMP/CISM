# Local pFUnit files.
set(pf_sources First_tests.pf)

# Sources to test.
set(sources_needed
  glad_constants.F90 glimmer_global.F90 glimmer_physcon.F90)
extract_sources("${sources_needed}" "${cism_sources}" test_sources)

# Do source preprocessing and add the executable.
add_pfunit_ctest(First_tests
  TEST_SOURCES "${pf_sources}"
  OTHER_SOURCES "${test_sources}")
  # Can rename pf_sources and test_sources to whatever I want later






