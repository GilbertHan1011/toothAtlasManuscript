Annotation
===============================

In this section, we will adapt the annotation pipeline from `scHarmonization <https://github.com/lsteuernagel/scHarmonization/>`_ :cite:p:`steuernagelHypoMapUnifiedSinglecell2022` to our data.
You can just run this pipeline directly. Here, I just break down the pipeline into several steps for better understanding.

.. toctree::
    :caption: annotation
    :maxdepth: 2
    :glob:

    20241025_mes_cluster
    step1_choose_optimal_cluster_level
    step2_mtree
    step3_marker_detection
    step4_prune_tree
    step5_annotation



Details of Each Level Annotation
--------------------------------

In this section, we explain the details of each level annotation. 

Annotate tooth atlas is a difficult task, bacause of following reasons:

- There is no standardized or authoritative annotation system for dental cell types.
- Annotation markers often conflict across different studies, making it difficult to establish consistent cell type definitions.
- Many existing annotation markers are specific to particular developmental stages or spatial locations, limiting their utility in integrated atlas studies.
- Current annotations predominantly rely on spatial concepts rather than molecular signatures, which may not accurately reflect distinct cell types and their biological functions.

Here are our basic annotation strategy:

- We ensure every annotation level is exhaustive.
- On this basis, we strive to use marker genes that are consistent across different studies.

You can explore the marker genes of each cluster in the `shiny app <https://zyflab.shinyapps.io/tooth/>`_.

.. toctree::
    :caption: Cluster details of mesenchyme
    :maxdepth: 2

    Level1_anno
    Level2_anno
    Level3_anno
    Level4_anno
    Level5_anno

.. toctree::
    :caption: Cluster details of epithelium
    :maxdepth: 2

    Level1_anno_epi
    Level2_anno_epi
    Level3_anno_epi