import sys

# Set a custom exception for user input errors to avoid showing traceback in the terminal
class InputError (Exception):
    pass

# Set a custom exception handler where our input error exception has a quiet behaviour
def custom_excepthook (exception, message, traceback):
    # Quite behaviour if it is our input error exception
    if exception == InputError:
        print('{0}: {1}'.format(exception.__name__, message))  # Only print Error Type and Message
        return
    # Default behaviour otherwise
    sys.__excepthook__(exception, message, traceback)
sys.excepthook = custom_excepthook