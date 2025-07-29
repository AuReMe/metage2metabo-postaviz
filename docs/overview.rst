Overview
========

M2M-PostAViz (_M2M Post-Analysis and Visualization_) is an interactive exploration platform for the exploration of metabolic potential predictions performed by [Metage2Metabo (M2M)](https://github.com/AuReMe/metage2metabo/tree/main).

M2M predicts reachable metabolites by communities of microorganisms by computing a Boolean abstraction (Network expansion) of metabolic activity. It takes as input a collection of genome-scale metabolic networks (GSMNs) and a description of the nutritional environment, then computes the molecules reachable by each GSMN individually, and by the community of GSMNs. In the latter case, it takes into account mutualistic interactions (cross-feeding) that may occur within the community, and therefore increase its overall metabolic potential.

Several outputs can be distinguished: what is produced by each member alone (_individual scope_), what becomes producible only through metabolic interactions (_added-value_), and the union of both (_community scope_). For each compound in one of these metabolic potentials, the information related to the producer(s) is retrieved by M2M and provided in dedicated tables. Therefore, we retrieve in some parts of the application the distinction between the individual metabolic potentials and the community metabolic potential, as well as the study of the cooperation added-value in some analyses. When nothin is mentioned, the application considers the community metabolic potential.

.. image:: ./pictures/m2m_overview.png
   :alt: General overview of Metage2Metabo's pipeline
   :width: 70%

M2M-PostAViz integrates and analyses all these outputs, especially in a context when multiple runs of M2M are performed, each one aiming at studying the metabolic potential of a microbial community. A typical use-case would be to run M2M for a cohort of microbiome samples, each described by a collection of GSMNs, for instance. In a cohort, the data comes with metadata that is used by M2M-PostAViz to analyse M2M's results and explore whether the predicted metabolic potential is statistically associated with metadata variables.

.. image:: ./pictures/postaviz_overview.png
   :alt: General overview of M2M-PostAViz
   :width: 70%

License
=======

GNU Lesser General Public License v3 (LGPLv3)

Authors
=======

Léonard Brindel and `Clémence Frioux <https://cfrioux.github.io>`__ -- `Inria Pleiade team <https://team.inria.fr/pleiade/>`__

Acknowledgements
================

- David James Sherman
- Jean-Marc Frigerio
- Pleiade team members
