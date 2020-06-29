import logging

logger = logging.getLogger('minimap')
if logger.handlers:
    logger.handlers = []

# Remove root log messages
# logger = logging.getLogger()
# original_logging_handlers = logger.handlers[:]
# for handler in original_logging_handlers:
#     if isinstance(handler, logging.FileHandler):
#         logger.removeHandler(handler)

logger.setLevel(logging.DEBUG)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter(
    '{levelname} {asctime} {name} {lineno}: {message}', style='{')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

logger.propagate = False