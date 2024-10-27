selector_to_html = {"a[href=\"#motivation\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Motivation<a class=\"headerlink\" href=\"#motivation\" title=\"Permalink to this heading\">#</a></h2><p>In level 2 annotation, we devided the mesenchyme into 9 clusters, including four large populations (\u201cC9-1\u201d,\u201dC9-2\u201d,\u201dC9-6\u201d,\u201dC9-7\u201d) and five small populations. In this level, we want to assign the consensus annotation from literature to these clusters.\n<img alt=\"png\" src=\"../_images/anno_C9.png\"/></p>", "a[href=\"#level-2-annotation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Level 2 annotation<a class=\"headerlink\" href=\"#level-2-annotation\" title=\"Permalink to this heading\">#</a></h1><h2>Motivation<a class=\"headerlink\" href=\"#motivation\" title=\"Permalink to this heading\">#</a></h2><p>In level 2 annotation, we devided the mesenchyme into 9 clusters, including four large populations (\u201cC9-1\u201d,\u201dC9-2\u201d,\u201dC9-6\u201d,\u201dC9-7\u201d) and five small populations. In this level, we want to assign the consensus annotation from literature to these clusters.\n<img alt=\"png\" src=\"../_images/anno_C9.png\"/></p>", "a[href=\"#c9-1-bud-mesenchyme\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">C9-1 : Bud Mesenchyme<a class=\"headerlink\" href=\"#c9-1-bud-mesenchyme\" title=\"Permalink to this heading\">#</a></h2><p>In the study focused on tooth germ development <span id=\"id1\">[<a class=\"reference internal\" href=\"../references.html#id27\" title=\"Yaofeng Wang, Yifan Zhao, Shubin Chen, Xiaoming Chen, Yanmei Zhang, Hong Chen, Yuansong Liao, Jiashu Zhang, Di Wu, Hongxing Chu, Hongying Huang, Caixia Wu, Shijuan Huang, Huichao Xu, Bei Jia, Jie Liu, Bo Feng, Zhonghan Li, Dajiang Qin, Duanqing Pei, and Jinglei Cai. Single cell atlas of developing mouse dental germs reveals populations of CD24+ and Plac8+ odontogenic cells. Science Bulletin, 67(11):1154\u20131169, June 2022. doi:10.1016/j.scib.2022.03.012.\">Wang <em>et al.</em>, 2022</a>]</span> <span id=\"id2\">[<a class=\"reference internal\" href=\"../references.html#id2\" title=\"Junjun Jing, Jifan Feng, Yuan Yuan, Tingwei Guo, Jie Lei, Fei Pei, Thach-Vu Ho, and Yang Chai. Spatiotemporal single-cell regulatory atlas reveals neural crest lineage diversification and cellular function during tooth morphogenesis. Nature Communications, 13(1):4803, August 2022. doi:10.1038/s41467-022-32490-y.\">Jing <em>et al.</em>, 2022</a>]</span>, the data reveals significant heterogeneity among the datasets at different stages, such as the bud stage and cap stage. We have also observed that the thickening stage mesenchyme (&lt; E13.5) and Bud stage (E13.5-E14.5) only exist in the C9-1 cluster. And C9-2 cluster mostly composited by mesenchyme at the cap stage (E14.5-E18.5). Therefore, we assigned the C9-1 cluster as \u201cBud Mesenchyme\u201d and C9-2 as \u201cCap Mesenchyme\u201d.</p>"}
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
