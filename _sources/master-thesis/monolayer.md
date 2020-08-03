---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.9'
    jupytext_version: 1.5.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Hoi boi

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from myst_nb import glue
```

```{code-cell} ipython3
:tags: [remove-output, hide-input]

x = np.linspace(0, 1, 100)
fig, ax  = plt.subplots()
ax.plot(x, x**2)
glue("my-fig", fig, display=False)
```

bla bla

```{glue:figure} my-fig
---
figclass: margin-caption
name: figfig
---
Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer accumsan tellus nec pulvinar lacinia. Pellentesque sed tristique augue. Aenean in congue lorem. Duis commodo auctor est quis pharetra. Aenean interdum at ligula id blandit. Vestibulum feugiat dignissim leo, eu maximus ex convallis sed. Sed dolor lectus, viverra non nulla et, hendrerit vulputate nisi. Etiam nec orci laoreet, facilisis nunc id, semper ex.
```
bla bla

```{figure} ../logo.png
---
figclass: margin-caption
name: figfig2
---
Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer accumsan tellus nec pulvinar lacinia. Pellentesque sed tristique augue. Aenean in congue lorem. Duis commodo auctor est quis pharetra. Aenean interdum at ligula id blandit. Vestibulum feugiat dignissim leo, eu maximus ex convallis sed. Sed dolor lectus, viverra non nulla et, hendrerit vulputate nisi. Etiam nec orci laoreet, facilisis nunc id, semper ex.
```

+++

Dit was de functie:

$$
f(x) = x^2
$$
