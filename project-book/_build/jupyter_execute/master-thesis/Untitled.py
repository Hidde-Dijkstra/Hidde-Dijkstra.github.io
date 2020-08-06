import numpy as np
import matplotlib.pyplot as plt
from myst_nb import glue

fig, ax = plt.subplots()
x = np.linspace(0, 1, 100)
ax.plot(x, x**2)
glue("square_fig", fig, display=False)

```{glue:figure} square_fig
---
figclass: margin-caption
name: square
---
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.
```

