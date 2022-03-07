# SDS
SDS forms microstructures (micelles) at most accessible simulation conditions. This presents some interesting questions in its modeling. 

This case study uses the following coarse-grained mapping:

![SDS](SDS.png){: style="height:89px;width:250px"}

where every 2 tail carbon beads are grouped to a single bead, and the entire $SO_4^-$ head group is mapped to a single bead.

The case study also builds off of [water](../water/index.md), [dodecane](../dodecane/index.md), and [NaCl](../salt/index.md) interactions derived following previous case studies. As outlined in the [theoretical considerations]("bestpractices.md"), this takes some care in choosing a coarse graining ensemble. For reproducibile structuring we deposit SDS at a dodecane-water interface:

![dodecane-water-SDS](dodecane-water-SDS.png){: style="height:106px;width:334px"}

(The bounding boxes above outline the dodecane region sandwiched between two SDS interfacial layers.)

See the [jupyter notebook](SDS.ipynb) implementing this in coarse graining.