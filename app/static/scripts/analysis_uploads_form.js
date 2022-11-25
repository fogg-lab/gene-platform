function cancelReqListener() {
    // file upload successfully cancelled
    console.log(this.responseText);
}

function updateNextButton() {
    // hide next button if progress counts or coldata not uploaded
    document.getElementById("next_button").disabled = !(isUploaded("counts") && isUploaded("coldata"));
}

function goto(url) {
    window.location.href = url;
}
