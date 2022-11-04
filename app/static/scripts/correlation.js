var UPLOAD_ENDPOINT = "/upload";
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
    cancel("counts.tsv");
    update_next_button();
}

function cancel(filename) {
    var cancelReq = new XMLHttpRequest();
    let filename_base = filename.split(".")[0];
    let progress_bar_id = "progress_of_" + filename_base;
    let upload_div_id = filename_base + "_upload_div";
    let cancel_button_id = "cancel_" + filename_base + "_button";
    let cancel_req_query = "filename=" + filename;
    let file_progress_bar = document.getElementById(progress_bar_id);
    cancelReq.addEventListener("load", cancelReqListener);
    cancelReq.open("POST", "/cancel-upload");
    cancelReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");

    cancelReq.send(cancel_req_query);

    // hide the progress bar
    file_progress_bar.remove();
    progress_div = document.getElementById(progress_bar_id + "_div");
    progress_div.style.backgroundPosition = "100% 0";
    progress_div.style.display = 'none';

    // hide the cancel button
    document.getElementById(cancel_button_id).style.display = 'none';

    // show the file upload div
    document.getElementById(upload_div_id).style.display = 'block';
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