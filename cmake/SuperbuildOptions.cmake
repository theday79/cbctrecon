##-----------------------------------------------------------------------------
##  See COPYRIGHT.TXT and LICENSE.TXT for copyright and license information
##-----------------------------------------------------------------------------
##  These macros are used to forward variables defined in the outer build
##  to the inner build.  The variable values cannot be a list; these
##  don't expand properly.
##-----------------------------------------------------------------------------
# macro: sb_variable
# Add a variable to the list of superbuild variables that are passed to the
# inner build
macro (sb_variable _var)
  list (APPEND sb_cmake_vars ${_var})
endmacro ()

# macro: sb_set
# Set a variable and add it to the list of superbuild variables
# that are passed to the inner build
macro (sb_set _var)
  set (${_var} ${ARGN})
  list (APPEND sb_cmake_vars ${_var})
endmacro ()

# macro: sb_option
# Create an option in the cmake-gui, and mark the variable as a superbuild
# variable to be passed to the inner build
macro (sb_option _var _desc _defval)
  option (${_var} ${_desc} ${_defval})
  list (APPEND sb_cmake_vars ${_var})
endmacro ()

# macro: sb_option_enum
# Create an enum option in the cmake-gui, and mark the variable as a
# superbuild variable to be passed to the inner build
macro (sb_option_enum _var _desc _defval)
  option_enum (${_var} ${_desc} ${_defval} ${ARGN})
  list (APPEND sb_cmake_vars ${_var})
endmacro ()
