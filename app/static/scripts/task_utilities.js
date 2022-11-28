var overlay_timer = {};
var fetch_task_updates = {};
var seconds_elapsed = 0;
var no_progress_count = 0;
var last_log_offset = 0;

function cancelTask() {
    task_cancelled = true;
    overlay_div = getOverlay();
    overlay_div.style.display = "none";
    cancel_task_button = document.getElementById("cancel_task_button")
    cancel_task_button.remove()
    clearInterval(overlay_timer);
    clearInterval(fetch_task_updates);
    overlay_timer = {};
    seconds_elapsed = 0;
    document.getElementById("time").innerHTML = "";
}


const escapeHtml = (unsafe) => {
    return unsafe.replaceAll('&', '').replaceAll('<', '').replaceAll('>', '').replaceAll('"', '&quot;').replaceAll("'", '&#039;');
}


function taskFailed(status_msg) {
    overlay = document.getElementById("overlay");
    overlay.remove();
    let messages = document.getElementById("messages");
    messages.innerHTML = status_msg;
}


function taskLogReqListener() {
    let log_with_metadata = this.responseText.split("###PROGRESS_METADATA:");
    let log_update_text = log_with_metadata[0].slice(0, -1);    // remove trailing newline
    // metadata is a JSON string
    let log_metadata = JSON.parse(log_with_metadata[1]);
    last_log_offset = log_metadata["last_log_offset"];
    let task_status = log_metadata["task_status"];

    if (task_status == "completed") {
        loadResults();
        return
    }

    if (log_update_text.length == 0) {
        no_progress_count += 1;
        let timeout = 120;
        let timeout_notice_interval = 30;
        if (no_progress_count % timeout_notice_interval == 0 && no_progress_count < timeout) {
            log_update_text = `No console output for ${no_progress_count} seconds. `;
            log_update_text += `Timeout=${timeout} seconds.\n`;
        }
        else if (no_progress_count == timeout) {
            taskFailed("Error: No response from server.");
            clearInterval(overlay_timer);
            clearInterval(fetch_task_updates);
            overlay_timer = {};
            seconds_elapsed = 0;
        }
        return;
    }

    no_progress_count = 0;

    task_log = document.getElementById("task_log");
    if (task_log != null) {
        let console_text = document.getElementById("console_text");
        console_text_scroll_pos_y = -1;
        console_text_scroll_pos_x = -1;
        if (console_text != null) {
            if (console_text.scrollTop < (console_text.scrollHeight - console_text.clientHeight) - 20) {
                console_text_scroll_pos_y = console_text.scrollTop;
            }
            console_text_scroll_pos_x = console_text.scrollLeft;
        }
        let task_log_html = task_log.innerHTML;
        let pre_open_tag = "<pre id='console_text' style='width: 100%; height: 100%'>"
        if ((task_log_html.length > pre_open_tag.length + 6) && task_log_html.slice(-6) == "</pre>") {
            task_log_html = `${pre_open_tag}${escapeHtml(task_log_html.slice(pre_open_tag.length,-6) + log_update_text)}</pre>`;
        } else {
            task_log_html = `${pre_open_tag}${escapeHtml(log_update_text)} </pre>`;
        }
        if (task_log_html.length > 18000) {
            new_start_index = task_log_html.length - 15000;
            task_log_html = `${pre_open_tag}${escapeHtml(task_log_html.slice(new_start_index).slice(0,-6))}</pre>`;
        }
        task_log.innerHTML = task_log_html;
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


function taskUpdateListener() {
    let progress_metadata = this.responseText.split("###PROGRESS_METADATA:")[1];
    let task_status = JSON.parse(progress_metadata)["task_status"];
    if (task_status == "completed") {
        loadResults();
    }
}


function loadResults() {
    clearInterval(overlay_timer);
    clearInterval(fetch_task_updates);
    overlay_timer = {};
    seconds_elapsed = 0;
    task_id = getTaskID();
    overlay = getOverlay();
    overlay.remove();
    displayResults();
}


function getTaskID() {
    let taskIDElement = document.getElementById("taskID");
    let taskID = "";
    if (taskIDElement) {
        taskID = taskIDElement.value;
    }
    return taskID;
}


function getTaskUpdate(log_offset=0) {
    var taskLogReq = new XMLHttpRequest();
    console_log = document.getElementById("console_log");
    if (console_log != null && console_log.style.display != "none") {
        taskLogReq.addEventListener("load", taskLogReqListener);
    }
    else {
        taskLogReq.addEventListener("load", taskUpdateListener);
    }
    task_id = getTaskID();
    dest = `get-progress?task_id=${task_id}&last_log_offset=${log_offset}`;
    taskLogReq.open("GET", dest);
    taskLogReq.send();
}


function startTaskUpdateRepeater() {
    fetch_task_updates = setInterval(function() {
        getTaskUpdate(last_log_offset);
    }, 1000);
}


function showRuntime() {
    clearInterval(overlay_timer);
    overlay_timer = {};
    seconds_elapsed = 0;
    overlay_div = getOverlay();
    console_log = document.getElementById("console_log");
    console_log.style.display="block";
    task_cancelled = false;
    const cancel_task_button = document.createElement("button");
    cancel_task_button.id = "cancel_task_button";
    cancel_task_button.innerHTML = "Cancel";
    cancel_task_button.onclick = cancelTask;
    overlay_div.appendChild(cancel_task_button);
    cancel_task_button.className = "button";
    overlay_timer = setInterval(function () {
        time_elem = document.getElementById("time")
        if (time_elem != null) {
            time_elem.innerHTML = "Time elapsed: " + seconds_elapsed + " seconds";
        }
        seconds_elapsed += 1;
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

    clearInterval(overlay_timer);
    clearInterval(fetch_task_updates);

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
                <div id="task_log"></div>
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
