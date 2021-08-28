import logging.config

DETAILED_FORMAT = "%(asctime)s %(levelname)s %(pathname)s:%(lineno)d: %(message)s"


def init_logger(logfile=None, level=logging.INFO, logger_name=None):
    """
    Initialize LOGGER and set logfile

    Args:
        logfile (str):
        level:
        logger_name (str):

    """

    if logger_name:
        logger = logging.getLogger(logger_name)
    else:
        logger = logging.getLogger()

    logger.setLevel(level)
    ch = logging.StreamHandler()
    if logfile:
        fh = logging.FileHandler(logfile)
    else:
        fh = None

    for handler in (x for x in (ch, fh) if x is not None):
        handler.setFormatter(logging.Formatter(DETAILED_FORMAT))
        logger.addHandler(handler)

    return logger
