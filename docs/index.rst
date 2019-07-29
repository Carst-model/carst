.. carst documentation master file, created by
   sphinx-quickstart on Mon Jul 22 15:52:22 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the documentation for the Carst Model!
=================================================

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   carst


See the firedrake documentation at https://firedrakeproject.org/firedrake.html for a complement to these documents.

Example Usage
=============

Using Diamond
-------------

Using diamond provides a user-friendly way to utilise carst. Diamond can be started with "make diamond" in the root directory of the model. Save the document produced as *diamond_input.xml* in the root dir. Carst will run automatically when diamond closes.

Using Python Scripts
--------------------

Using python scripts, the full functionality of *carst* can be utilised. The user should prepare the initial setting for the model:

* **The mesh the model will operate in**. This must be an object of type *firedrake.mesh.MeshGeometry*.

  A simple rectangular mesh could be initialised with:

  .. code-block:: python

      my_mesh = fd.RectangleMesh(50, 25, 10000, 5000)

* **The layout of the seabed**. This must be in the form of a callable taking two arguments:

  * The coordinate space of the mesh the model is working on. This will be passed to it as a firedrake object (indexable), of type *firedrake.ufl.geometry.SpacialCoordinate*.

  * The function space of the mesh the model is working on. This will be passed to it as a firedrake object, of type *firedrake.functionspaceimpl.WithGeometry*.

  The callable must return an object of type *firedrake.function.Function*; we recommend using:

  .. code-block:: python

     firedrake.project(<your function>, function_space)

  For example, to implement the following expression:

  :math:`100 \tanh(\frac{1}{2000} (x_0 - 6000))`

  We would compose the following python code:

  .. code-block:: python

     def my_land(coordinate_space, function_space):
         return fd.project(100 * fd.tanh(0.0005 * (coordinate_space[0] - 6000)), function_space)

  Note that any mathematical functions used are implemented in the form of the provided functions in *firedrake*. See the firedrake documentation for more details.

* **The function describing the sea level**. This must be passed as a string to be evaluated by python's **eval()**.

  An example:

  .. code-block:: python

     sea_level_lit = "25 * fd.sin(T / 50000 * 180 / math.pi)"

* **The times the model is to function on**. These are passed as numeric types in a 3-length tuple as follows, composed of:

  .. code-block:: python

     times = (<start time>, <time step>, <output time>)

* **The output folder**. This is a simple relative path to mark the folder which files should be output into, passed as a string.

  .. code-block:: python

      output_folder = "output"

* **Which processes should be enabled**. At the time of writing, only "diffusion" and "carbonates" are implemented. These should be bools.

  .. code-block:: python

      diff_enabled = True
      carbs_enabled = False

* **Extra options**. At the time of writing, "diffusion" requires a diffusion coefficient to be given, a float named "diff_coeff".

  .. code-block:: python

      my_diff_coeff = 2.0

Once all these options are gathered, they can be passed (in the above order) to a *carst.options.CarstOptions* constructor with the *carst.options.initialisation_method.raw_values* ini_type:

.. code-block:: python

   from carst import options, solver

   my_options = options.CarstOptions(
       options.initialisation_method.raw_values,
       my_mesh,
       my_land,
       sea_level_lit,
       times,
       output_folder,
       diffusion=diff_enabled,
       carbonates=carbs_enabled,
       diff_coeff=my_diff_coeff,
   )

The *CarstOptions* class should be passed to the constructor of *carst.solver.CarstModel*:

.. code-block:: python

   my_solver = solver.CarstModel(my_options)

Finally, an initial condition for the sediment must be specified. This, similarly to the seabed function, must be passed as a callable which returns an object of type *firedrake.function.Function*.

In order to implement:

:math:`\frac{10000}{\sqrt{2 \pi \times 250^2}} e^{- \frac{(x_0 - 6000)^2}{2 \times 250^2}} + \frac{25000}{\sqrt{2 \pi \times 1000^2}} e^{- \frac{(x_0 - 4000)^2}{2 \times 1000^2}}`

We would use:

.. code-block:: python

   def my_initial_cond(coordinate_space, function_space):
       return fd.project(
           (
               10000 / (fd.sqrt(2 * math.pi * 250**2))
               * fd.exp(-((coordinate_space[0] - 6000)**2) / (2 * 250**2))
           ) + (
               25000 / (fd.sqrt(2 * math.pi * 1000**2))
               * fd.exp(-((coordinate_space[0] - 4000)**2) / (2 * 1000**2))
           ),
           function_space
       )

This must then be passed to the *set_condition* method of your *CarstModel* instance:

.. code-block:: python

   my_solver.set_condition(my_initial_cond)

Finally, the model is ready to begin simulation. *CarstModel.advance()* should be used to advance the model by 1 time step at a time.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
