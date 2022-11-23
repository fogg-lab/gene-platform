// hide progress bar divs and cancel buttonsfor files which have not been uploaded
if (document.getElementById("progress_of_counts") == null) {
    document.getElementById("progress_of_counts_div").style.display = 'none';
    document.getElementById("cancel_counts_button").style.display = 'none';
} else {
    document.getElementById("progress_of_counts_div").style.backgroundPosition = "0% 0";
    document.getElementById("counts_upload_div").style.display = 'none';
}

// disable next button if progress counts not uploaded
update_next_button()

function cancelReqListener() {
    // file upload successfully cancelled
    console.log(this.responseText);
}

function cancel_counts() {
    cancel_upload("counts.tsv");
    update_next_button();
}

function update_next_button() {
    // hide next button if progress counts or coldata not uploaded
    if (
        document.getElementById("progress_of_counts") == null
        || document.getElementById("progress_of_counts").className != "success"
        ) {
            document.getElementById("next_button").disabled = true;
        }
        else {
            document.getElementById("next_button").disabled = false;
        }
}

function corrReqListener() {
    // Done computing correlations
    let messages = document.getElementById("messages");
    messages.innerHTML = this.responseText;
    getPlots();
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

function submitRNASeqCorrelation() {
    corr_select = document.getElementById("corr_method")
    let corr_method = corr_select.options[corr_select.selectedIndex].value;
    let submitCorrReq = new XMLHttpRequest();
    let submit_corr_query = "corr_method=" + corr_method
    submitCorrReq.addEventListener("load", corrReqListener);
    submitCorrReq.open("POST", "/submit-rnaseq-correlation");
    submitCorrReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    submitCorrReq.send(submit_corr_query);
    let messages = document.getElementById("messages");
    messages.innerHTML = "Computing sample correlations...";
}

function getPlots() {
    for (corr_method of ["pearson", "spearman"]) {
        let getPlotsReq = new XMLHttpRequest();
        let getPlotsQuery = "?corr_method=" + corr_method;
        if (corr_method == "pearson") {
            getPlotsReq.addEventListener("load", getPearsonPlotListener);
        } else {
            getPlotsReq.addEventListener("load", getSpearmanPlotListener);
        }
        dest = `/get-correlation-plot${getPlotsQuery}`;
        getPlotsReq.open("GET", dest);
        getPlotsReq.send();
    }
}