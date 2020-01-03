import logging
import sys

def setup_logger(level='INFO',logfile=None):

    root = logging.getLogger()

    screen_log_format = '[%(levelname).4s] [%(funcName)25s] %(message)s'
    file_log_format = '[%(asctime)-15s] [%(levelname)08s]  [%(name)s] [%(funcName)s] %(message)s'

    if level.upper() not in ['DEBUG','INFO','WARNING','ERROR','CRITICAL']:
        level_value=logging.WARNING        
    if level == 'DEBUG':
        level_value = logging.DEBUG
    if level == 'INFO':
        level_value = logging.INFO
    if level == 'WARNING':
        level_value = logging.WARNING
    if level == 'ERROR':
        level_value = logging.ERROR
    if level == 'CRITICAL':
        level_value = logging.CRITICAL

    root.handlers = []

    screen_handler = logging.StreamHandler(sys.stdout)
    screen_handler.setLevel(level_value)
    screen_handler.setFormatter(logging.Formatter(screen_log_format))

    root.addHandler(screen_handler)

    if logfile is not None:
        print("Logging to file not implemented yet.")

    return()
    
