################################################################################
# Current Date
################################################################################

# For the moment, we can't use string(timestamp...) for getting the current date,
# because it needs a 2.8.11 version for cmake. Since we may have to compile the 
# stuff on servers where we are not admin, we may face too old cmake versions for
# supporting this feature. 
# Ideally, we should use: string (TIMESTAMP gatb-core-date "%Y-%m-%d %H:%M:%S")

# The trick is to rely on PERL (likely to be installed on the machine)
# (see http://osdir.com/ml/programming.tools.cmake.user/2006-06/msg00323.html)

INCLUDE(FindPerl)

MACRO (GetCurrentDate dateString)
    # We execute a command that retrieves the current date.
    EXECUTE_PROCESS (
        COMMAND "${PERL_EXECUTABLE}" "-le" "@T=localtime; printf (\"%04d-%02d-%02d/%02d:%02d:%02d\",$T[5]+1900,$T[4]+1,$T[3],$T[2],$T[1],$T[0])"
        OUTPUT_VARIABLE ${dateString}
    )
ENDMACRO()


MACRO (GetCurrentDateShort dateString)
    # We execute a command that retrieves the current date.
    EXECUTE_PROCESS (
        COMMAND "${PERL_EXECUTABLE}" "-le" "@T=localtime; printf (\"%04d%02d%02d\",$T[5]+1900,$T[4]+1,$T[3])"
        OUTPUT_VARIABLE ${dateString}
    )
ENDMACRO()
