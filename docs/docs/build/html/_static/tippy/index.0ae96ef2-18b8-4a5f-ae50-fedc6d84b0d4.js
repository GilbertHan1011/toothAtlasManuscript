selector_to_html = {"a[href=\"trajectory/index.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Trajectory<a class=\"headerlink\" href=\"#trajectory\" title=\"Permalink to this heading\">#</a></h1><p>In this section, we focus on odontogenesis and amelogenesis -\ntwo critical developmental processes in tooth formation.\nOur aim is to explore and understand the dynamic patterns of gene expression\nduring these processes. A key question we seek to address is whether there exists a\nuniversal molecular axis that governs mineralization across different contexts.</p>", "a[href=\"#welcome-to-scatlas-s-documentation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Welcome to scAtlas\u2019s documentation!<a class=\"headerlink\" href=\"#welcome-to-scatlas-s-documentation\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>ScAtlas is a pipeline for building large-scale single-cell atlases. While there is exponential growth\nin the amount of single-cell data generated, a high quality reference atlas can be valuable resource for\nresearch groups to derive insights from their data. Here, we hope to provide a pipeline for building\nlarge-scale single-cell atlases that can be used for a variety of applications.</p>", "a[href=\"annotation/index.html#details-of-each-level-annotation\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Details of Each Level Annotation<a class=\"headerlink\" href=\"#details-of-each-level-annotation\" title=\"Permalink to this heading\">#</a></h2><p>In this section, we explain the details of each level annotation.</p><p>Annotate tooth atlas is a difficult task, bacause of following reasons:</p>", "a[href=\"#contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Contents<a class=\"headerlink\" href=\"#contents\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"trajectory/multipotent.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Stemness and Start Point Inference<a class=\"headerlink\" href=\"#stemness-and-start-point-inference\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>Determining the starting point of a trajectory presents a significant challenge, especially when working with integrated datasets. To address this, we can leverage key biological indicators such as cellular age and developmental potential. In this section, we will employ multiple computational methods including Cytotrace2 <span id=\"id1\">[<a class=\"reference internal\" href=\"references.html#id39\">Kang <em>et al.</em>, 2024</a>]</span>, Cytotrace <span id=\"id2\">[<a class=\"reference internal\" href=\"references.html#id43\" title=\"Gunsagar S. Gulati, Shaheen S. Sikandar, Daniel J. Wesche, Anoop Manjunath, Anjan Bharadwaj, Mark J. Berger, Francisco Ilagan, Angera H. Kuo, Robert W. Hsieh, Shang Cai, Maider Zabala, Ferenc A. Scheeren, Neethan A. Lobo, Dalong Qian, Feiqiao B. Yu, Frederick M. Dirbas, Michael F. Clarke, and Aaron M. Newman. Single-cell transcriptional diversity is a hallmark of developmental potential. Science (New York, N.Y.), 367(6476):405, January 2020. doi:10.1126/science.aax0249.\">Gulati <em>et al.</em>, 2020</a>]</span>, and SCENT <span id=\"id3\">[<a class=\"reference internal\" href=\"references.html#id44\" title=\"Andrew E. Teschendorff and Tariq Enver. Single-cell entropy for accurate estimation of differentiation potency from a cell's transcriptome. Nature Communications, 8(1):15599, June 2017. doi:10.1038/ncomms15599.\">Teschendorff and Enver, 2017</a>]</span> to systematically infer the developmental potential of cells and identify suitable trajectory starting points.</p>", "a[href=\"annotation/index.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Annotation<a class=\"headerlink\" href=\"#annotation\" title=\"Permalink to this heading\">#</a></h1><p>In this section, we will adapt the annotation pipeline from <a class=\"reference external\" href=\"https://github.com/lsteuernagel/scHarmonization/\">scHarmonization</a> <span id=\"id1\">[<a class=\"reference internal\" href=\"references.html#id22\" title=\"Lukas Steuernagel, Brian Y. H. Lam, Paul Klemm, Georgina K. C. Dowsett, Corinna A. Bauder, John A. Tadross, Tamara Sotelo Hitschfeld, Almudena del Rio Martin, Weiyi Chen, Alain J. de Solis, Henning Fenselau, Peter Davidsen, Irene Cimino, Sara N. Kohnke, Debra Rimmington, Anthony P. Coll, Andreas Beyer, Giles S. H. Yeo, and Jens C. Br\u00fcning. HypoMap\u2014a unified single-cell gene expression atlas of the murine hypothalamus. Nature Metabolism, 4(10):1402\u20131419, October 2022. doi:10.1038/s42255-022-00657-y.\">Steuernagel <em>et al.</em>, 2022</a>]</span> to our data.\nYou can just run this pipeline directly. Here, I just break down the pipeline into several steps for better understanding.</p>", "a[href=\"#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>ScAtlas is a pipeline for building large-scale single-cell atlases. While there is exponential growth\nin the amount of single-cell data generated, a high quality reference atlas can be valuable resource for\nresearch groups to derive insights from their data. Here, we hope to provide a pipeline for building\nlarge-scale single-cell atlases that can be used for a variety of applications.</p>"}
skip_classes = ["headerlink", "sd-stretched-link"]

window.onload = function () {
    for (const [select, tip_html] of Object.entries(selector_to_html)) {
        const links = document.querySelectorAll(`div.content ${select}`);
        for (const link of links) {
            if (skip_classes.some(c => link.classList.contains(c))) {
                continue;
            }

            tippy(link, {
                content: tip_html,
                allowHTML: true,
                arrow: true,
                placement: 'auto-start', maxWidth: 500, interactive: false,
                onShow(instance) {MathJax.typesetPromise([instance.popper]).then(() => {});},
            });
        };
    };
    console.log("tippy tips loaded!");
};