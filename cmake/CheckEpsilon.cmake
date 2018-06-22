##---------------------------------------------------------------------------
## See COPYRIGHT.TXT and LICENSE.TXT for copyright and license information
##---------------------------------------------------------------------------

macro (check_epsilon _out_var)
  try_run (RUN_RESULT_VAR COMPILE_RESULT_VAR
    ${PLM_BINARY_DIR}
    ${PLM_SOURCE_DIR}/cmake/test_eps.cxx
    RUN_OUTPUT_VARIABLE ${_out_var})
endmacro ()
