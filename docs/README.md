to build: `rm -rf _build && make html`
to open: `xdg-open _build/html/index.html`

for math in docstrings, see [this](https://documentation.help/Sphinx/math.html)
Basically, you need to double `\\` the latex specials, or use python (`r""`) raw strings

