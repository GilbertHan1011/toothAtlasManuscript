selector_to_html = {"a[href=\"#introduction\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>ScAtlas is a pipeline for building large-scale single-cell atlases. While there is exponential growth\nin the amount of single-cell data generated, a high quality reference atlas can be valuable resource for\nresearch groups to derive insights from their data. Here, we hope to provide a pipeline for building\nlarge-scale single-cell atlases that can be used for a variety of applications.</p>", "a[href=\"pre-integration/scib-benchmark.html#why-this-matters\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Why this matters<a class=\"headerlink\" href=\"#why-this-matters\" title=\"Permalink to this heading\">#</a></h2><p>For atlas level scRNA analysis, integration algrithoms are essential to keep biological signals while removing technical noises <span id=\"id1\">[<a class=\"reference internal\" href=\"references.html#id19\" title=\"Malte D. Luecken, M. B\u00fcttner, K. Chaichoompu, A. Danese, M. Interlandi, M. F. Mueller, D. C. Strobl, L. Zappia, M. Dugas, M. Colom\u00e9-Tatch\u00e9, and Fabian J. Theis. Benchmarking atlas-level data integration in single-cell genomics. Nature Methods, 19(1):41\u201350, January 2022. doi:10.1038/s41592-021-01336-8.\">Luecken <em>et al.</em>, 2022</a>]</span>. The performance of integration methods highly depends on the data itself <span id=\"id2\">[<a class=\"reference internal\" href=\"references.html#id19\" title=\"Malte D. Luecken, M. B\u00fcttner, K. Chaichoompu, A. Danese, M. Interlandi, M. F. Mueller, D. C. Strobl, L. Zappia, M. Dugas, M. Colom\u00e9-Tatch\u00e9, and Fabian J. Theis. Benchmarking atlas-level data integration in single-cell genomics. Nature Methods, 19(1):41\u201350, January 2022. doi:10.1038/s41592-021-01336-8.\">Luecken <em>et al.</em>, 2022</a>]</span>. Therefore we recommend to benchmark integration methods before you start integration.</p><p>Here, we use <a class=\"reference external\" href=\"https://github.com/theislab/scib-pipeline\">scib benchmarking pipeline</a> to perform integration benchmark.</p>", "a[href=\"#welcome-to-scatlas-s-documentation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Welcome to scAtlas\u2019s documentation!<a class=\"headerlink\" href=\"#welcome-to-scatlas-s-documentation\" title=\"Permalink to this heading\">#</a></h1><h2>Introduction<a class=\"headerlink\" href=\"#introduction\" title=\"Permalink to this heading\">#</a></h2><p>ScAtlas is a pipeline for building large-scale single-cell atlases. While there is exponential growth\nin the amount of single-cell data generated, a high quality reference atlas can be valuable resource for\nresearch groups to derive insights from their data. Here, we hope to provide a pipeline for building\nlarge-scale single-cell atlases that can be used for a variety of applications.</p>", "a[href=\"#contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Contents<a class=\"headerlink\" href=\"#contents\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"pre-integration/scib-benchmark.html#process\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Process<a class=\"headerlink\" href=\"#process\" title=\"Permalink to this heading\">#</a></h2><p><a class=\"reference external\" href=\"https://github.com/theislab/scib-pipeline\">scib-pipeline</a> offers a handy snakemake pipeline to perform integration benchmark. However, it lacks maintainance since two years ago. You might have to take some times to fix the pacakge dependencies conflicts.</p><p>The process is as follows:</p>", "a[href=\"pre-integration/scib-benchmark.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Integration methods benchmarking<a class=\"headerlink\" href=\"#integration-methods-benchmarking\" title=\"Permalink to this heading\">#</a></h1><h2>Why this matters<a class=\"headerlink\" href=\"#why-this-matters\" title=\"Permalink to this heading\">#</a></h2><p>For atlas level scRNA analysis, integration algrithoms are essential to keep biological signals while removing technical noises <span id=\"id1\">[<a class=\"reference internal\" href=\"references.html#id19\" title=\"Malte D. Luecken, M. B\u00fcttner, K. Chaichoompu, A. Danese, M. Interlandi, M. F. Mueller, D. C. Strobl, L. Zappia, M. Dugas, M. Colom\u00e9-Tatch\u00e9, and Fabian J. Theis. Benchmarking atlas-level data integration in single-cell genomics. Nature Methods, 19(1):41\u201350, January 2022. doi:10.1038/s41592-021-01336-8.\">Luecken <em>et al.</em>, 2022</a>]</span>. The performance of integration methods highly depends on the data itself <span id=\"id2\">[<a class=\"reference internal\" href=\"references.html#id19\" title=\"Malte D. Luecken, M. B\u00fcttner, K. Chaichoompu, A. Danese, M. Interlandi, M. F. Mueller, D. C. Strobl, L. Zappia, M. Dugas, M. Colom\u00e9-Tatch\u00e9, and Fabian J. Theis. Benchmarking atlas-level data integration in single-cell genomics. Nature Methods, 19(1):41\u201350, January 2022. doi:10.1038/s41592-021-01336-8.\">Luecken <em>et al.</em>, 2022</a>]</span>. Therefore we recommend to benchmark integration methods before you start integration.</p><p>Here, we use <a class=\"reference external\" href=\"https://github.com/theislab/scib-pipeline\">scib benchmarking pipeline</a> to perform integration benchmark.</p>", "a[href=\"pre-integration/scib-benchmark.html#results\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Results<a class=\"headerlink\" href=\"#results\" title=\"Permalink to this heading\">#</a></h2>", "a[href=\"#indices-and-tables\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Indices and tables<a class=\"headerlink\" href=\"#indices-and-tables\" title=\"Permalink to this heading\">#</a></h1>"}
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