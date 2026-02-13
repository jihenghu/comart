"""
Definition of the Exception specific to COMART.
"""

import warnings


class ComartError(Exception):
    """
    Error raised by the model.
    """

    pass


class ComartWarning(Warning):
    """
    Warning raised by the model.
    """

    pass


def comart_warn(message, **kwargs):
    """
    Issue a COMART warning with a standard message about disabling warnings.
    """

    disable_warning_message = """To disable all comart warnings, use: import warnings;
warnings.filterwarnings("ignore", comae.error.ComartWarning). See the warnings Python documentation for finer
controls.

"""
    if "stacklevel" in kwargs:
        kwargs["stacklevel"] += 1
    else:
        kwargs["stacklevel"] = 2
    warnings.warn(message + "\n\n" + disable_warning_message, category=ComartWarning, **kwargs)