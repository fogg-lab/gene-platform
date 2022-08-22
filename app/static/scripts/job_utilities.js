var timer = {};
var counter = 0;
no_log_update_count = 0;


function cancelJob() {
    job_cancelled = true;
    overlay_div = getOverlay();
    overlay_div.style.display = "none";
    cancel_job_button = document.getElementById("cancel_job_button")
    cancel_job_button.remove()
    clearInterval(timer);
    timer = {};
    counter = 0;
    document.getElementById("time").innerHTML = "";
}


const escapeHtml = (unsafe) => {
    return unsafe.replaceAll('&', '').replaceAll('<', '').replaceAll('>', '').replaceAll('"', '&quot;').replaceAll("'", '&#039;');
}


function jobFailed(status_msg) {
    overlay = document.getElementById("overlay");
    overlay.remove();
    let messages = document.getElementById("messages");
    messages.innerHTML = status_msg;
}


function jobLogReqListener() {
    log_update_text = this.responseText;
    if (log_update_text.length == 0) {
        no_log_update_count += 1;
        let timeout = 120;
        let timeout_notice_interval = 30;
        if (no_log_update_count % timeout_notice_interval == 0 && no_log_update_count < timeout) {
            log_update_text = `No console output for ${no_log_update_count} seconds. `;
            log_update_text += `Timeout=${timeout} seconds.\n`;
        }
        else if (no_log_update_count == timeout) {
            jobFailed("Error: No response from server.");
            clearInterval(timer);
            timer = {};
            counter = 0;
            return;
        }
    }
    else {
        no_log_update_count = 0;
    }
    job_log = document.getElementById("job_log");
    if (job_log != null) {
        let console_text = document.getElementById("console_text");
        console_text_scroll_pos_y = -1;
        console_text_scroll_pos_x = -1;
        if (console_text != null) {
            if (console_text.scrollTop < (console_text.scrollHeight - console_text.clientHeight) - 20) {
                console_text_scroll_pos_y = console_text.scrollTop;
            }
            console_text_scroll_pos_x = console_text.scrollLeft;
        }
        let job_log_html = job_log.innerHTML;
        let pre_open_tag = "<pre id='console_text' style='width: 100%; height: 100%'>"
        if ((job_log_html.length > pre_open_tag.length + 6) && job_log_html.slice(-6) == "</pre>") {
            job_log_html = `${pre_open_tag}${escapeHtml(job_log_html.slice(pre_open_tag.length,-6) + log_update_text)}</pre>`;
        } else {
            job_log_html = `${pre_open_tag}${escapeHtml(log_update_text)} </pre>`;
        }
        if (job_log_html.length > 18000) {
            new_start_index = job_log_html.length - 15000;
            job_log_html = `${pre_open_tag}${escapeHtml(job_log_html.slice(new_start_index).slice(0,-6))}</pre>`;
        }
        job_log.innerHTML = job_log_html;
        console_text = document.getElementById("console_text");
        if (console_text_scroll_pos_y != -1) {
            console_text.scrollTop = console_text_scroll_pos_y;
        } else {
            console_text.scrollTop = console_text.scrollHeight;
        }
        if (console_text_scroll_pos_x != -1) {
            console_text.scrollLeft = console_text_scroll_pos_x;
        }
    }
}


function requestConsoleOutput() {
    overlay = document.getElementById("overlay");
    if (overlay != null && overlay.style.display != "none") {
        var jobLogReq = new XMLHttpRequest();
        jobLogReq.addEventListener("load", jobLogReqListener);
        jobLogReq.open("GET", "/get-console-output");
        jobLogReq.send();
    }
}


function showRuntime() {
    clearInterval(timer);
    timer = {};
    counter = 0;
    overlay_div = getOverlay();
    console_log = document.getElementById("console_log");
    console_log.style.display="block";
    job_cancelled = false;
    const cancel_job_button = document.createElement("button");
    cancel_job_button.id = "cancel_job_button";
    cancel_job_button.innerHTML = "Cancel";
    cancel_job_button.onclick = cancelJob;
    overlay_div.appendChild(cancel_job_button);
    cancel_job_button.className = "button";
    timer = setInterval(function () {
        time_elem = document.getElementById("time")
        if (time_elem != null) {
            time_elem.innerHTML = "Time elapsed: " + counter + " seconds";
        }
        counter += 1;
        requestConsoleOutput();
    }, 1000);
}


function confirmSubmission(confirmation_text) {
    dialog_box = document.getElementById("confirm_div_id")
    if (dialog_box == null) {
        dialog_box = document.createElement("div");
    }
    dialog_box.id = "confirm_div_id";
    dialog_box.innerHTML = confirmation_text;
    dialog_box.innerHTML += "<button id='okay_btn_id'>Okay</button>\n";
    dialog_box.innerHTML += "\n<button id='cancel_btn_id'>Cancel</button>\n";
    dialog_box.style.display="block";

    clearInterval(timer);

    overlay = getOverlay();
    
    console_log = document.getElementById("console_log");
    console_log.style.display="none";
    overlay.appendChild(dialog_box)
    document.getElementById('okay_btn_id').onclick = function() {
        dialog_box.style.display="none";
        overlay.style.display="none";
        if (!(confirmation_text.includes("Error"))) {
            submissionConfirmed();
        }
    };
    document.getElementById('cancel_btn_id').onclick = function() {
        dialog_box.style.display="none";
        overlay.style.display="none";
    };
}


function getOverlay() {
    let overlay = document.getElementById("overlay");
    if (overlay == null) {
        overlay = document.createElement("div");
        overlay.id = "overlay";
        overlay.innerHTML = `
            <span id="time"></span>
            <details id="console_log">
                <summary style="font-size: 20px; font-weight:bold;">Show console output</summary>
                <div id="job_log"></div>
            </details>
        `;
        document.body.appendChild(overlay);
    }
    overlay.style.display="flex";
    return overlay;
}


function submitReqListener() {
    var confirmation_text = this.responseText;
    confirmSubmission(confirmation_text);
}
