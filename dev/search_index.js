var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "#CatmullRom.jl-1",
    "page": "Overview",
    "title": "CatmullRom.jl",
    "category": "section",
    "text": ""
},

{
    "location": "#Centripetal-parameterization-for-Catmull-Rom-interpoint-connectivity.-1",
    "page": "Overview",
    "title": "Centripetal parameterization for Catmull-Rom interpoint connectivity.",
    "category": "section",
    "text": ""
},

{
    "location": "#Copyright-and-thinsp;2018-2019-by-Jeffrey-Sarnoff.-and-nbsp;-and-nbsp;-This-work-is-released-under-The-MIT-License.-1",
    "page": "Overview",
    "title": "Copyright ©&thinsp;2018-2019 by Jeffrey Sarnoff. &nbsp;&nbsp;  This work is released under The MIT License.",
    "category": "section",
    "text": ""
},

{
    "location": "#[![Build-Status](https://travis-ci.org/JeffreySarnoff/CatmullRom.jl.svg?branchmaster)](https://travis-ci.org/JeffreySarnoff/CatmullRom.jl)-1",
    "page": "Overview",
    "title": "(Image: Build Status)",
    "category": "section",
    "text": ""
},

{
    "location": "perspective/#",
    "page": "Perspecive",
    "title": "Perspecive",
    "category": "page",
    "text": ""
},

{
    "location": "perspective/#Perspective-1",
    "page": "Perspecive",
    "title": "Perspective",
    "category": "section",
    "text": "Catmull-Rom splines are a workhorse of computer graphics. Using the centripetal parameterization, they become a very handy general purpose tool for fast, attractive curvilinear blending. Often, they give interpoint \"motion\" a naturalistic feel.(Image: hair)<p align=\"center\"><a href=\"http://www.cemyuksel.com/research/catmullrom_param\">Cem Yuksel\'s Research</a></p><p align=\"center\"><a href=\"http://www.cemyuksel.com/research/catmullromparam/catmullromcad.pdf\">Parameterization and Applications of Catmull-Rom Curves</a></p>(Image: CR_Centripetal)   (Image: OneCurve123)<p align=\"center\">The blue curves show Catmull-Rom paths that obtain using the centripetal parameterization (α=0.5)</p>"
},

{
    "location": "threefunctions/#",
    "page": "Three Functions",
    "title": "Three Functions",
    "category": "page",
    "text": "Centripetal Catmull-Rom Pathmaking\n\ncatmullrom( points )\n\ncatmullrom_by_arclength( points )\n"
},

{
    "location": "threefunctions/#Uniform-Intermediation-1",
    "page": "Three Functions",
    "title": "Uniform Intermediation",
    "category": "section",
    "text": "catmullrom( pointsalonga_path )\ncatmullrom( pointsalongapath, narcsbetweenneighbors )There are two ways to connect path-adjacent points using Centripetal Catmull-Rom machinery. The most often used iterplaces a given number of curvilinear waypoints between each adjacent pair of original points.  All neighbors become connected by that given number of intermediating places. Though the places differ, the proportional advancing between abcissae is consistent. There is a default for this intermediation count, nonetheless trying a few different values may help you visualize the significance that is of import.crpoints = catmullrom( points )\n\ncrpoints = catmullrom( points, n_arcs_between_neighbors )"
},

{
    "location": "threefunctions/#Arclength-Relative-Allocation-1",
    "page": "Three Functions",
    "title": "Arclength Relative Allocation",
    "category": "section",
    "text": "When the points\' coordinates are spread differently along distinct axes, the interpoint distances along one coordinate have a very different nature from the intercoordinate spreads along another coordinate axis.  The distances separating adjacent point pairs may vary substantively.  This is particularly true when working in higher dimensional regions of an orthonormal coordinate space.  One may use more intermediating placements between adjacent points that are relatively far apart, and fewer between adjacent points that are in close relative proximity.crpoints = catmullrom_byarc( points )\n\ncrpoints = catmullrom_byarc( points, (min_arcs_between_points, max_arcs_between_points) )using CatmullRom, Plots\n\nresult = catmullrom(points, n_arcs_per_pair)    # your points, how many new points to place between adjacents\n                                                 # result is a vector of coordinates, e.g. [xs, ys, zs]\nplot(result...,)When your points have nonuniform separation, or separation extents vary with coordinate dimension, it is of benefit to allocate more of the new inbetween points where there are relatively greater distances between your adjacent points.  The most appropriate measure for comparison and weighting is interpoint arclength.  This package implements a well-behaved approximation to Catmull-Rom arclengths appropriate to the centripetal parameterization.  You can use this directly.using CatmullRom, Plots\n\nresult = catmullrom_by_arclength(points) # result is a vector of coordinates, e.g. [xs, ys, zs]\n \nresult = catmullrom_by_arclength(points, (atleast_min_arcs_total, atmost_max_arcs_total))\n                                          # min, max pertain to the whole path of the curve\nxs, ys = result\nplot(xs, ys)"
},

{
    "location": "threefunctions/#Open-and-Closed-Curves-1",
    "page": "Three Functions",
    "title": "Open and Closed Curves",
    "category": "section",
    "text": "CatmullRom processes the extremal points of closed curves differently from open curves. A curve in which the first and last points are identical is recognized as closed. A function, close_seq!, is available to ensure curves intended to be closed are made closed in an exact and proper way. It is good practice to use this function with closed curves, and so assure they are crisp where the escribed path rejoins itself.close_seq!( points )            # this is the only function that may change some part of your data\n                                # any change is limited to copying the first point into the last \npoints = close_seq!( points )   # (the same thing)"
},

{
    "location": "Examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "Examples/#Centripetal-Catmull-Rom-Examples-sup[𝓪](#source)/sup-1",
    "page": "Examples",
    "title": "Centripetal Catmull-Rom Examples <sup>𝓪</sup>",
    "category": "section",
    "text": "|                     shape                               |              detail                           |   |:––––––––––––––––––––––––––––:|:–––––––––––––––––––––––-:|   |                                                          |                                                 |   | <img src=\"https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/CatmullRomcircledpihalf.png\" width=\"300\">  |      <img src=\"https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/CatmullRomsectionofcircle.png\" width=\"300\">|   |                                                          |                                                 |    |    <img src=\"https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/teardropbyarc.png\" width=\"400\">          |   <img src=\"https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/teardropsection.PNG\" width=\"250\">         |    |    <img src=\"https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/sphericalspiral.png\" width=\"400\">          |   <img src=\"https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/sphericalspirallocus.png\" width=\"130\">    |  <p align=\"center\"><sup><a name=\"source\">𝓪</sup> <a href=\"https://github.com/JeffreySarnoff/CatmullRom.jl/tree/master/examples\">using this package to generate some of these examples</a></p>  "
},

{
    "location": "pointsalongapath/#",
    "page": "Points along a path",
    "title": "Points along a path",
    "category": "page",
    "text": ""
},

{
    "location": "pointsalongapath/#Points-along-a-path-1",
    "page": "Points along a path",
    "title": "Points along a path",
    "category": "section",
    "text": "A sequence of 2D, 3D .. nD points is required.  There is no limit on the number of coordinate dimensions.   The first coordinate of each point become the abcissae (e.g. the x coordinate values).  The second [, third etc.] become [successive] ordinates (e.g. the ys, zs ...).Every point in a givne sequence must has the same number of constiuent coordinates.  Coordinates are considered to be values along orthonormal axes.  All ordinate axes are fitted with respect to the same abcissae. So, the arcs that connect successive ys are arcs hewn from a succession of (x_i, y_i) ordered pairs and the arcs connecting successive zs are arcs hewn from a succession of (x_i, z_i) ordered pairs.  It is easy to work with other axial pairings. To generate arcs using the sequence of (y_i, z_i) pairs: ys_zs = catmullrom( collect(zip(ys, zs)) ).The point sequence itself may be provided as a vector of points or as a tuple of points. Type used for a Point example coordinates are retrievable you support\n   \nsmall vector [1.0, 3.5 ] coord(point, i) = point[i] builtin\nsmall tuple (1.0, 3.5) coord(point, i) = point[i] builtin\n   \nStaticVector SVector( 1.0, 3.5 ) coord(point, i) = point[i] builtin\nNamedTuple (x = 1.0, y = 3.5 ) coord(point, i) = point[i] builtin\n   \nstruct Point(1.0, 3.5) coord(point, i) = point[i] getindex\n   struct Point{T}\n    x::T\n    y::T\n    z::T\nend\n\nfunction Base.getindex(point::Point{T}, i::Integer) where T\n    if i == 1\n       point.x\n    elseif i == 2\n       point.y\n    elseif i == 3\n       point.z\n    else\n       throw(DomainError(\"i must be 1, 2, or 3 (not $i)\"))\n    end\nend"
},

{
    "location": "hints/#",
    "page": "A few hints",
    "title": "A few hints",
    "category": "page",
    "text": ""
},

{
    "location": "hints/#hints-1",
    "page": "A few hints",
    "title": "hints",
    "category": "section",
    "text": "If your points are disaggregated (e.g. all the xs in vecofxs, all the ys in vecofys)\naggregate them this way points = collect(zip(xs, ys, zs))\nOften, the xs respect x[i-1] < x[i] < x[i+1] or x[i-1] > x[i] > x[i+1]     \nwhen the path is a closed curve, one triplet may follow    \nx[i-1] < x[i] > x[i+1] or x[i-1] > x[i] < x[i+1]"
},

{
    "location": "references/#",
    "page": "References",
    "title": "References",
    "category": "page",
    "text": ""
},

{
    "location": "references/#references-1",
    "page": "References",
    "title": "references",
    "category": "section",
    "text": "Parameterization and Applications of Catmull-Rom CurvesThe Centripetal Catmull-Rom SplineCatmull-Rom spline without cusps or self-intersections"
},

]}
