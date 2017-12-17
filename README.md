This repo is forked from https://github.com/ISA-tools/isa-api/tree/py2_isatools-lite and needs some updating.

<img align="left" src="http://www.isa-tools.org/wp-content/uploads/2016/03/isa-api-logo.png" width="150px">
<br>

----
<img align="right" src="http://www.isa-tools.org/wp-content/themes/isatools-wp-theme/img/footer_logo.svg">
The open source ISA metadata tracking tools help to manage an increasingly diverse set of life science, environmental and biomedical experiments that employing one or a combination of technologies.

Built around the ‘Investigation’ (the project context), Study’ (a unit of research) and ‘Assay’ (analytical measurement) general-purpose Tabular format, the ISA tools helps you to provide rich description of the experimental metadata (i.e. sample characteristics, technology and measurement types, sample-to-data relationships) so that the resulting data and discoveries are reproducible and reusable.

To find out more about ISA, see [www.isa-tools.org](http://www.isa-tools.org)

To find out who's using ISA and about the ISA development and user community, see [www.isacommons.org](http://www.isacommons.org)

The *ISA API*  aims to provide you, the developer, with a set of tools to help you easily and quickly build your own ISA objects, validate, and convert between serializations of ISA-formatted datasets and other formats/schemas (e.g. SRA schemas). The ISA API is published on PyPI as the `isatools` package.

This project is the `isatools-rwval` package that provides ISA format read/write/validation functionality that the ISA API builds on.

Thie project contains the following modules:

 - `isatools.model` contains classes implementing the ISA Abstract Model as Python objects.
 - `isatools.isatab` contains features for parsing and serializing the ISA-Tab format to/from Python objects. The module also contains a ISA-Tab validator.
 - `isatools.isajson` contains features for parsing and serializing the ISA-JSON format to/from Python objects. The module also contains a ISA-JSON validator.
 - `isatools.isaviz` contains features for rendering ISA Python objects in visual artifacts with `matplotlib`.

 These modules should be considered the gold-standard utilities for using the ISA formats with the Python programming language.

[![Build Status](https://travis-ci.org/ISA-tools/isa-rwval.svg?branch=master)](https://travis-ci.org/ISA-tools/isa-rwval)
[![Coverage Status](https://coveralls.io/repos/github/ISA-tools/isa-rwval/badge.svg?branch=master)](https://coveralls.io/github/ISA-tools/isa-rwval?branch=master)
[![Documentation Status](https://readthedocs.org/projects/isatools/badge/?version=latest)](http://isatools.readthedocs.org/en/latest/?badge=latest)

----
*Authors*: [Code contributors](https://github.com/ISA-tools/isatools-core/graphs/contributors).

*License*: This code is licensed under the [? License](https://raw.githubusercontent.com/ISA-tools/isatools-core/master/LICENSE.txt).

*Repository*: [https://github.com/ISA-tools/isatools-core](https://github.com/ISA-tools/isatools-core)

*ISA team email*: [isatools@googlegroups.com](mailto:isatools@googlegroups.com)

*ISA discussion group*: [https://groups.google.com/forum/#!forum/isaforum](https://groups.google.com/forum/#!forum/isaforum)

*Github issue tracker*: [https://github.com/ISA-tools/isa-api/issues](https://github.com/ISA-tools/isatools-core/issues)
