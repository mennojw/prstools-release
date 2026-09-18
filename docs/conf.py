

import json
from configparser import ConfigParser
from pathlib import Path

from sphinxawesome_theme.postprocess import Icons


ROOT = Path(__file__).resolve().parents[1]

settings = ConfigParser()
settings.read(ROOT / "settings.ini")
metadata = settings["DEFAULT"]

project = metadata.get("title", "prstools")
author = metadata.get("author", "Menno Witteveen et al.")
release = metadata.get("version", "")
version = release


extensions = [
    "myst_nb",
    "sphinx_togglebutton",
        "sphinxarg.ext",
     "sphinx_copybutton",
]

source_suffix = {
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
}

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**/.ipynb_checkpoints",
]
suppress_warnings = ["misc.highlighting_failure"]
templates_path = ["_templates"]
html_static_path = ["_static"]
html_css_files = ["custom.css"]

# html_theme = "sphinxawesome_theme"
html_theme = "sphinx_rtd_theme"
html_title = f"prstools {release} documentation" if release else "prstools documentation"
html_theme_options = {
    "main_nav_links": {
        "GitHub": "https://github.com/mennojw/prstools-release",
        "PyPI": "https://pypi.org/project/prstools/",
    },
    "awesome_external_links": True,
    "awesome_headerlinks": True,
    "show_prev_next": True,
}
html_permalinks_icon = Icons.permalinks_icon

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]

# Render the outputs already stored in the notebooks. This keeps local and
# Read the Docs builds fast and avoids requiring tutorial data during a build.

smartquotes = False

nb_execution_mode = "off"
nb_merge_streams = True

nb_code_prompt_show = "Show command output"
nb_code_prompt_hide = "Hide command output"

# nb_execution_mode = "off"
# nb_merge_streams = True
# nb_code_prompt_show = "OUTPUT — Show full output"
# nb_code_prompt_hide = "OUTPUT — Hide output"


# Automatically collapse substantial shell-command output in notebooks while
# leaving the command itself visible. Notebook authors do not need cell tags.
cli_output_collapse_after = 100

# 
def _content_size(value):
    if isinstance(value, str):
        return len(value)
    if isinstance(value, dict):
        return sum(_content_size(item) for item in value.values())
    if isinstance(value, list):
        return sum(_content_size(item) for item in value)
    return 0


def _collapse_cli_outputs(app, docname, source):
    notebook_path = Path(app.srcdir, f"{docname}.ipynb")
    if not notebook_path.is_file():
        return

    notebook = json.loads(source[0])
    changed = False
    for cell in notebook.get("cells", []):
        command = "".join(cell.get("source", [])).lstrip()
        outputs = cell.get("outputs", [])
        if (
            cell.get("cell_type") == "code"
            and command.startswith("!")
            and _content_size(outputs) >= cli_output_collapse_after
        ):
            tags = cell.setdefault("metadata", {}).setdefault("tags", [])
            if "hide-output" not in tags:
                tags.append("hide-output")
                changed = True

    if changed:
        source[0] = json.dumps(notebook)


def setup(app):
    app.connect("source-read", _collapse_cli_outputs)
