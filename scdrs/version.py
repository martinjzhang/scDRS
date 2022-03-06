# Store the version here so:
# 1) we don't load dependencies by storing it in __init__.py
# 2) we can import it in setup.py for the same reason
# 3) we can import it into your module module
__version__ = '1.0.1'
__version_info__ = tuple([ int(num) for num in __version__.split('.')])