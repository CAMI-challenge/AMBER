import os
import os.path as pt
from scripttest import TestFileEnvironment

root_dir = pt.abspath(pt.join(pt.dirname(__file__), '..'))

def before_scenario(context, _):
    context.env = TestFileEnvironment(base_path = context.config.userdata["TMPDIR"])

    # Mounting volumes in Docker needs an explict full path.
    # The path cannot be hard coded into the features as it varies
    # Therefore dynamically pass the full path prefix as the TMPDIR
    # environment variable.
    for (key, value) in context.config.userdata.items():
        os.environ[key] = value
