var UPLOAD_ENDPOINT = "/upload";

// hide progress bar divs and cancel buttons for files which have not been uploaded
if (document.getElementById("progress_of_counts") == null) {
    document.getElementById("progress_of_counts_div").style.display = 'none';
    document.getElementById("cancel_counts_button").style.display = 'none';
} else {
    document.getElementById("progress_of_counts_div").style.backgroundPosition = "0% 0";
    document.getElementById("counts_upload_div").style.display = 'none';
}

if (document.getElementById("progress_of_coldata") == null) {
    document.getElementById("progress_of_coldata_div").style.display = 'none';
    document.getElementById("cancel_coldata_button").style.display = 'none';
} else {
    document.getElementById("progress_of_coldata_div").style.backgroundPosition = "0% 0";
    document.getElementById("coldata_upload_div").style.display = 'none';
}

// disable next button if progress counts or coldata not uploaded
update_next_button()

function cancelReqListener() {
    // file upload successfully cancelled
    console.log(this.responseText);
}

function bcReqListener() {
    // batch correction completed
    // display download button
    let messages = document.getElementById("messages");
    messages.innerHTML = this.responseText;
    if (this.responseText.includes("Batch correction complete.")) {
        document.getElementById("download-bc-btn").disabled = false;
    }
}

function cancel_counts() {
    cancel("counts.tsv");
    update_next_button();
}

function cancel_coldata() {
    cancel("coldata.tsv");
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
        || document.getElementById("progress_of_coldata") == null
        || document.getElementById("progress_of_counts").className != "success"
        || document.getElementById("progress_of_coldata").className != "success"
        ) {
            document.getElementById("next_button").disabled = true;
        }
        else {
            document.getElementById("next_button").disabled = false;
        }
}

function submit_batch_correction() {
    datatype_select = document.getElementById("data_type")
    let data_type = datatype_select.options[datatype_select.selectedIndex].text
    let reference_level = document.getElementById("reference_level").value;
    let contrast_level = document.getElementById("contrast_level").value;
    let submitbcReq = new XMLHttpRequest();
    let submit_bc_query = "data_type=" + data_type + "&reference_level=" + reference_level + "&contrast_level=" + contrast_level;
    submitbcReq.addEventListener("load", bcReqListener);
    submitbcReq.open("POST", "/submit-batch-correction");
    submitbcReq.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    submitbcReq.send(submit_bc_query);
    let messages = document.getElementById("messages");
    document.getElementById("download-bc-btn").disabled = true;
    messages.innerHTML = "Performing batch correction...";
}
