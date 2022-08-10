from pathlib import Path
import inspect
import sys
from matplotlib.testing.compare import compare_images
import matplotlib.testing.exceptions


def caller_name(skip=2):
    """Get a name of a caller in the format module.class.method

    `skip` specifies how many levels of stack to skip while getting caller
    name. skip=1 means "who calls me", skip=2 "who calls my caller" etc.

    An empty string is returned if skipped levels exceed stack height
    """

    def stack_(frame):
        framelist = []
        while frame:
            framelist.append(frame)
            frame = frame.f_back
        return framelist

    stack = stack_(sys._getframe(1))
    start = 0 + skip
    if len(stack) < start + 1:  # pragma: no cover
        return ""  # pragma: no cover
    parentframe = stack[start]

    name = []
    module = inspect.getmodule(parentframe)
    # `modname` can be None when frame is executed directly in console
    # TODO(techtonik): consider using __main__
    if module:  # pragma: no branch
        name.append(module.__name__)
    # detect classname
    if "self" in parentframe.f_locals:
        # I don't know any way to detect call from the object method
        # XXX: there seems to be no way to detect static method call - it will
        #      be just a function call
        name.append(parentframe.f_locals["self"].__class__.__name__)  # pragma: no cover
    codename = parentframe.f_code.co_name
    if codename != "<module>":  # pragma: no branch  # top level usually
        name.append(codename)  # function or a method
    del parentframe
    return ".".join(name)


def caller_file(skip=2):
    """Get a name of a caller in the format module.class.method

    `skip` specifies how many levels of stack to skip while getting caller
    name. skip=1 means "who calls me", skip=2 "who calls my caller" etc.

    An empty string is returned if skipped levels exceed stack height
    """

    def stack_(frame):
        framelist = []
        while frame:
            framelist.append(frame)
            frame = frame.f_back
        return framelist

    stack = stack_(sys._getframe(1))
    start = 0 + skip
    if len(stack) < start + 1:  # pragma: no cover
        return ""  # pragma: no cover
    parentframe = stack[start]

    return parentframe.f_code.co_filename


def dump_cp_for_changed_images(generated_image_path, should_path):
    import shlex

    print("use %s to accept all image changes" % test_accept_image_path)
    if not test_accept_image_path.exists():
        test_accept_image_path.write_text("#!/usr/bin/sh\n")
    with open(test_accept_image_path, "a") as op:
        op.write(
            "cp %s %s\n"
            % (shlex.quote(str(generated_image_path)), shlex.quote(str(should_path)))
        )
    test_accept_image_path.chmod(0o755)


def assert_image_equal(generated_image_path, suffix="", tolerance=2, should_path=None):
    """assert that the generated image and the base_images/{test_module}/{cls}/{function_name}{suffix}.extension is identical"""
    generated_image_path = Path(generated_image_path).absolute()
    if should_path is None:
        extension = generated_image_path.suffix
        caller = caller_name(1)
        caller_fn = caller_file(1)
        parts = caller.split(".")
        if len(parts) >= 3:
            func = parts[-1]
            cls = parts[-2]
            module = parts[-3]
            # if cls.lower() == cls:  # not actually a class, a module instead
            # module = cls
            # cls = "_"
        else:
            module = parts[-2]
            cls = "_"
            func = parts[-1]
        should_path = (
            Path(caller_fn).parent
            / "base_images"
            / module
            / cls
            / (func + suffix + extension)
        ).resolve()
    if not generated_image_path.exists():
        raise IOError(f"Image {generated_image_path} was not created")
    if not should_path.exists():
        should_path.parent.mkdir(exist_ok=True, parents=True)
        dump_cp_for_changed_images(generated_image_path, should_path)
        raise ValueError(
            f"Base_line image not found, perhaps: \ncp {generated_image_path} {should_path}"
        )
    try:

        err = compare_images(
            str(should_path), str(generated_image_path), tolerance, in_decorator=False
        )
    except matplotlib.testing.exceptions.ImageComparisonFailure as e:
        dump_cp_for_changed_images(generated_image_path, should_path)
        raise ValueError(
            "Matplot lib testing for \n%s \n%s failed\n%s"
            % (generated_image_path, should_path, e)
        )
    # if isinstance(err, ValueError):
    # raise ValueError(
    # "Images differed significantly, rms: %.2f\nExpected: %s\n, Actual: %s\n, diff: %s\n Accept with cp %s %s\n"
    # % (err.rms, err.expected, err.actual, err.diff, err.actual, err.expected)
    # )
    if err is not None:
        dump_cp_for_changed_images(generated_image_path, should_path)
        raise ValueError(err)


if Path("run").exists():
    test_accept_image_path = Path("run/accept_all_image_changes.sh").absolute()
else:
    test_accept_image_path = Path("tests/run/accept_all_image_changes.sh").absolute()
if test_accept_image_path.exists():
    test_accept_image_path.unlink()
