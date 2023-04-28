---
# SPDX-FileCopyrightText: 2022-2023 CERN and the Allpix Squared authors
# SPDX-License-Identifier: CC-BY-4.0
title: "Neutron detection with Timepix3"
description: "This example demonstrate the possibility to simulate neutron detection by conversion using neutron capture media. A Timepix3 detector with a Boron10/Polyethylene thin coating is exposed to a source of neutrons with a energy of 0.025 eV / 10 MeV."
---

The physics list QGSP_BERT_HP is used. To obtain more realistic clusters for high energy deposit, the repulsion model of the `ProjectionPropagation` module, by setting the `repulsion_deposit` parameter as in the examples. The empirical `repulsion_attenuation_factor` allow to modify the intensity of the repulsion simulated. It is an empirical parameter to be tuned to your experimental situation. 


## Thermal neutrons with B10

Neutrons are captured by the Boron and a Lithium ion plus an alpha particle are emmited. The  detector will see one of the emmited particles, allowing detection.


## Fast neutron detection with polyethylene

For the polyethylene target, the neutron interact in the PE layer and a recoil proton is emmited and detected in the Silicon detector. 