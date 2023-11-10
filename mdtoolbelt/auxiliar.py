import sys

# Set a custom exception for user input errors to avoid showing traceback in the terminal
class QuietException (Exception):
    pass

# Set custom exception which are to be quiet because they are not caused by an error in the code
class InputError (QuietException):
    pass

class TestFailure (QuietException):
    pass

class MissingDependency (QuietException):
    pass

# Set a custom exception handler where our input error exception has a quiet behaviour
def custom_excepthook (exception_class, message, traceback):
    # Quite behaviour if it is our input error exception
    if QuietException in exception_class.__bases__:
        print('{0}: {1}'.format(exception_class.__name__, message))  # Only print Error Type and Message
        return
    # Default behaviour otherwise
    sys.__excepthook__(exception_class, message, traceback)
sys.excepthook = custom_excepthook