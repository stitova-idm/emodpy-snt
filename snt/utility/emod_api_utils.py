import emod_api.schema_to_class as s2c


def suppress_warnings(show_warnings: bool = True):
    """
    Turn on/off emod_api warnings.
    Args:
        show_warnings: True/False
    Returns:
        None
    """
    s2c.show_warnings = show_warnings
