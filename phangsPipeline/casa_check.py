def is_casa_installed():
    """Check if CASA is installed."""

    casa_enabled = False
    try:
        import casatasks
        import casatools
        casa_enabled = True
    except (ImportError, ModuleNotFoundError):
        pass

    return casa_enabled
