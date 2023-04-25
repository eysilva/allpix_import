# Ultra-Cold Neutron (UCN) detection Example

This example demonstrate the possibility to simulate neutron detection by conversion using neutron capture media. A Timepix3 detector with a Boron10 thin coating is exposed to a source of neutrons with a energy of 1e-9 eV. A second Timepix3 is placed on top of the first one, sensors facing back to back. T

he physics list QGSP_BERT_HP is used, in addition to the UCN special physics, activated by the `use_UCN_physics` flag.

Neutrons are captured by the Boron and a Lithium ion plus an alpha particle are emmited. Each detector detects one of the emmited particles. The Measurement of both decay products allow for a precise determination of the neutron origin. 

A Galodinium mask covers half of the detector. The sharpness of the transition in the neutron hitmap informs us on the neutron detection resolution of the assembly. 
