var task_cols = ["id", "type", "status", "created", "updated", "results", "remove"];
var next_sort_reverse = Object.fromEntries(task_cols.map( x => [x, false]));
var current_sorted_col = "updated";

function sort_tasks(th_id) {
    if (["results", "remove"].includes(th_id)) {
        return;
    }

    // get all tasks into an array of js objects
    var tasks_array = [];
    var table = document.getElementById("tasks-table");
    for (var i = 1; i < table.rows.length; i++) { // skip first row
        task = {
            "id":      table.rows[i].cells[0].innerHTML,
            "type":    table.rows[i].cells[1].innerHTML,
            "status":  table.rows[i].cells[2].innerHTML,
            "created": table.rows[i].cells[3].innerHTML,
            "updated": table.rows[i].cells[4].innerHTML,
            "results": table.rows[i].cells[5].innerHTML,
            "remove": table.rows[i].cells[6].innerHTML
        }
        tasks_array.push(task);
    }
    tasks_array.sort(function(a, b) {
        if (next_sort_reverse[th_id]) {
            return (b[th_id] > (a[th_id]) ?  1: -1);
        } else {
            return (a[th_id] > (b[th_id]) ?  1: -1);
        }
    })

    for (var i = 0; i < table.rows.length-1; i++) {
        for (var j = 0; j < table.rows[i].cells.length; j++) {
            field = task_cols[j];
            table.rows[i+1].cells[j].innerHTML = tasks_array[i][field];
        }
    }
    prev_sorted_th = document.getElementById(current_sorted_col);
    prev_sorted_th.classList.remove("sort-down");
    prev_sorted_th.classList.remove("sort-up");
    current_sorted_col = th_id;
    cur_sorted_th = document.getElementById(th_id);
    if (next_sort_reverse[th_id]) {
        cur_sorted_th.classList.add("sort-up");
    }
    else {
        cur_sorted_th.classList.add("sort-down");
    }
    next_sort_reverse[th_id] = !next_sort_reverse[th_id];
}

window.onload = function() {
    for (id of task_cols) {
        var dom_node = document.getElementById(id);
        dom_node.addEventListener("click", function() {
            sort_tasks(event.currentTarget.id);
        });
    }
    document.getElementById(current_sorted_col).classList.add("sort-up");
}
