var task_cols = ["id", "type", "status", "created", "updated", "results", "remove"];
var reverse = false

function sort_tasks(method) {
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
        if (reverse) {
            return (b[method] > (a[method]) ?  1: -1);
        } else {
            return (a[method] > (b[method]) ?  1: -1);
        }
    })

    for (var i = 0; i < table.rows.length-1; i++) {
        for (var j = 0; j < table.rows[i].cells.length; j++) {
            field = task_cols[j];
            table.rows[i+1].cells[j].innerHTML = tasks_array[i][field];
        }
    }
}

window.onload = function() {
    for (id of task_cols) {
        var dom_node = document.getElementById(id);
        dom_node.addEventListener("click", function() {
            sort_tasks(event.currentTarget.id);
            reverse = !reverse;
        })
    }
}
