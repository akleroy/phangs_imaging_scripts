
# Is CASA installed?

def is_casa_installed():

    casa_enabled = False

    # CASA 5
    try:
        import taskinit
        casa_enabled = True
        return casa_enabled
    except (ImportError, ModuleNotFoundError):
        pass

    # CASA 6
    try:
        import casatasks
        import casatools
        casa_enabled = True
        return casa_enabled
    except (ImportError, ModuleNotFoundError):
        pass

    return casa_enabled
