// hide progress bar divs and cancel buttonsfor files which have not been uploaded
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

if (document.getElementById("progress_of_filter") == null) {
    document.getElementById("progress_of_filter_div").style.display = 'none';
    document.getElementById("cancel_filter_button").style.display = 'none';
} else {document.getElementById("progress_of_filter_div").style.backgroundPosition = "0% 0";
    document.getElementById("filter_upload_div").style.display = 'none';
}

if (document.getElementById("progress_of_config") == null) {
    document.getElementById("progress_of_config_div").style.display = 'none';
    document.getElementById("cancel_config_button").style.display = 'none';
} else {document.getElementById("progress_of_config_div").style.backgroundPosition = "0% 0";
    document.getElementById("config_upload_div").style.display = 'none';
}

// disable next button if progress counts or coldata not uploaded
update_next_button()

function goto(url) {
    window.location.href = url;
}

function cancelReqListener() {
    // file upload successfully cancelled
    console.log(this.responseText);
}

function cancel_counts() {
    cancel_upload("counts.tsv");
    update_next_button();
}

function cancel_coldata() {
    cancel_upload("coldata.tsv");
    update_next_button();
}

function cancel_filter() {
    cancel_upload("filter.txt");
}

function cancel_config() {
    cancel_upload("config.yml");
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
