.. _alltasks:

List of tasks and parameters
############################

BBarolo's main algorithm for fitting 3D kinematic models to emission line data (:ref:`3DFIT <3dfit>`) makes use of a number of utilities. These tasks include, for example, the disk modeling (:ref:`GALMOD <galmodtask>`), the source finder (:ref:`SEARCH <searchtask>`) and the smoothing utility (:ref:`SMOOTH <smoothtask>`), and can be conveniently used outside the main algorithm as well.

In this page, I list the main tasks and related input parameters available in BBarolo. Parameter names are in **boldface**, default values are in brackets. The names of parameters are not case-sensitive.

.. toctree::
   :maxdepth: 1
   
   general
   fit3d
   spacepar
   galmod
   search
   smooth
   fit2d
   ellprof
   moment
