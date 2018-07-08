##---------------------------------------------------------------------------
## See COPYRIGHT.TXT and LICENSE.TXT for copyright and license information
##---------------------------------------------------------------------------

macro (check_char_sign _out_var)
  try_run (RUN_RESULT_VAR COMPILE_RESULT_VAR
    ${PLM_BINARY_DIR}
    ${PLM_SOURCE_DIR}/cmake/char_is_signed.cxx
    RUN_OUTPUT_VARIABLE ${_out_var})
endmacro ()
