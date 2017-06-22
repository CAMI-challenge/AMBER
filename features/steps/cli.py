import os
from behave import *

@when(u'I run the command')
def step_impl(context):
    context.output = context.env.run(
        "bash -c '{}'".format(os.path.expandvars(context.text)),
        expect_error  = True,
        expect_stderr = True)

def download_file(link, out):
    import wget
    return wget.download(link, out)

def get_stream(context, stream):
    assert stream in ['stderr', 'stdout'], "Unknown output stream {}".format(stream)
    return getattr(context.output, stream)

def assert_file_exists(file_):
    assert os.path.isfile(file_), "The file \"{}\" does not exist.".format(file_)

def get_env_path(context, file_):
    return os.path.join(context.env.cwd, file_)

def get_data_file_path(file_):
    dir_ = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(dir_, '..', '..', 'test', file_)

@then(u'the file "{}" should exist')
def step_impl(context, file_):
    assert_file_exists(get_env_path(context, file_))

@then(u'the exit code should be {code}')
def step_impl(context, code):
    returned = context.output.returncode
    assert returned == int(code), \
        "Process should return exit code {} but was {}".format(code, returned)

@given(u'I download the file "{link}" to "{dest}"')
def step_impl(context, link, dest):
    import sys
    normalized_dest = get_env_path(context, dest)
    sys.stdout = sys.__stdout__
    download_file(link, normalized_dest)

@given(u'I copy the example data files')
def step_impl(context):
    import shutil
    for row in context.table.rows:
        shutil.copy(get_data_file_path(row['source']),
            get_env_path(context, row['dest']))

@then(u'the {stream} should contain')
def step_impl(context, stream):
    output = get_stream(context, stream)
    assert context.text in output

@given(u'I downloaded the scripts')
def create_tmp_dir(context):
    dir_ = os.path.dirname(os.path.abspath(__file__))
    tmp_dir = os.path.join(dir_, '..', '..', "tmp")
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)