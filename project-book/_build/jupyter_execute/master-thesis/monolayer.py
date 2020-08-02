# Hoi boi

import numpy as np
import matplotlib.pyplot as plt
from myst_nb import glue

x = np.linspace(0, 1, 100)
fig, ax  = plt.subplots()
ax.plot(x, x**2)
glue("my-fig", fig, display=False)

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

Dit was de functie:

$$
f(x) = x^2
$$