# SIRF-SIMIND Connection

STIR user group meeting code that is slowly morphing into a more complete package. \

This will eventually be merged with https://github.com/samdporter/STIR-SIMIND-Connection \

I'll ensure better documentation shortly. In the meantime please see the demonstration notebook and python script (also callable by bash). I have included the Lutetium data on this repo. Please let me know if you would like access to the Y90 or Tc99m - I'll need to ask for permisson.

SIMIND can be downloaded [here](https://simind.blogg.lu.se/downloads/)
# Summary
This repo wraps SIMIND so that it can be configured with python. This functionality can be used on its own to avoid using the clunky (sorry, Michael) `change` fucntionality or can be used in tandem with the `SimindSimulator` class to simulate SIMIND projection data from SIRF `ImageData` and then output SIRF `AcquisitonData`

There is some code in its very early stages to use SIMIND as a forward projector within the SIRF framework. This should be used with caution, if at all. I'll get round to improving it soon.

# Tests
Currently there are no unit tests (I know, sorry). I and others have done some verification of SIMIND/STIR reconstructions using measured data and SIMIND itself has been investigated thouroughly. There's loads of literature out there, if you're interested.

# Note:
- .smc and .atn files must be in the same folder (atm input)

# ToDo:
- Update to use non-circular orbits using .cor file as template for SIRF
- Update to allow different kVp for attenuation corrections
- Loooots of documentation (READMEs + comments + docstrings)
