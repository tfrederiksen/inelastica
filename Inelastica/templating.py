from jinja2 import Environment, FileSystemLoader


j2 = Environment(loader=FileSystemLoader(["./", "/"]))


def _add_fdffloat_filter():
    def fdffloat(x, fmt="{:.8f}"):
        return fmt.format(x)
    j2.filters["fdffloat"] = fdffloat


_add_fdffloat_filter()
