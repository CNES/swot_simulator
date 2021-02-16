# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Exception handlers
------------------
"""
import linecache
import os
import traceback


def structured_traceback(exc: Exception,
                         call_stack: traceback.StackSummary) -> str:
    """Return a nice text describing the call stack which threw an
    exception.

    Args:
        exc (Exception): Exception raised.
        call_stack (traceback.StackSummary): Exception call stack thrown out.

    Returns:
        str: The text representing the call stack of the thrown exception.
    """
    exc_name = "%s.%s" % (exc.__module__, exc.__class__.__name__) if hasattr(
        exc, '__module__') else exc.__class__.__name__
    message = ["%s - Traceback (most recent call last):" % exc_name]

    for source, line_number, function, _ in call_stack:
        message.append("%s in %s" % (os.path.abspath(source), function))

        # Reading the Python code
        lines = linecache.getlines(source)

        # We will display the Python source code around the call so that
        # the user visualizes the incriminated lines well.
        start = line_number - 3
        end = line_number + 2
        start = max(start, 0)

        # Add the Python source code
        for idx, item in enumerate(lines[start:end]):
            where = start + idx + 1
            message.append(("----> " if where == line_number else "      ") +
                           str(where) + " " + item.rstrip())
        message.append("")

    # Finally, displays the error message
    message.append(exc_name + ": " + str(exc))
    return "\n".join(message)
