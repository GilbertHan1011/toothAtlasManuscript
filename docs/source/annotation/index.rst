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



Details of Each Level Annotation
--------------------------------

In this section, we explain the details of each level annotation. 
You can explore the marker genes of each cluster in the `shiny app <https://zyflab.shinyapps.io/tooth/>`_.

.. toctree::
    :caption: Cluster details
    :maxdepth: 2

    Level1_anno
    Level2_anno
    Level3_anno
    Level4_anno
    Level5_anno