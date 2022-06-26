function corrReqListener() {
    // Done computing correlations
    let messages = document.getElementById("messages");
    messages.innerHTML = this.responseText;
    get_plots();
}

function getPearsonPlotListener() {
    // Done getting plots
    let plot_div = document.getElementById("pearson_plot");
    plot_div.innerHTML = '<img src="data:image/png;base64,' + this.responseText + '"/>';
}

function getSpearmanPlotListener() {
    // Done getting plots
    let plot_div = document.getElementById("spearman_plot");
    plot_div.innerHTML = '<img src="data:image/png;base64,' + this.responseText + '"/>';
}

function submit_rnaseq_correlation() {
    corr_select = document.getElementById("corr_method")
    let corr_method = corr_select.options[corr_select.selectedIndex].value;
    let submitCorrReq = new XMLHttpRequest();
    let submit_corr_query = "corr_method=" + corr_method
    submitCorrReq.addEventListener("load", corrReqListener);
    submitCorrReq.open("POST", "/submit_rnaseq_sample_correlation");
    submitCorrReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    submitCorrReq.send(submit_corr_query);
    let messages = document.getElementById("messages");
    messages.innerHTML = "Computing sample correlations...";
}

function get_plots() {
    for (corr_method of ["pearson", "spearman"]) {
        let getPlotsReq = new XMLHttpRequest();
        let getPlotsQuery = "corr_method=" + corr_method;
        if (corr_method == "pearson") {
            getPlotsReq.addEventListener("load", getPearsonPlotListener);
        } else {
            getPlotsReq.addEventListener("load", getSpearmanPlotListener);
        }
        getPlotsReq.open("POST", "/get_" + corr_method + "_plot");
        getPlotsReq.send(getPlotsQuery);
    }
}