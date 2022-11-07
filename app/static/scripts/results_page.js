
display_filters = {
    "sortIndex": "",
    "search": "",
    "filter": false
}

displayed_data = [];
display_page_num = 1;

filter_selection = document.getElementById("filter_select");
sort_selection = document.getElementById("sort_select");
search_field = document.getElementById("search_field");


function populate_rows(start_index, end_index) {
    data = displayed_data.slice(start_index, end_index+1)
    new_output_tbody = document.createElement("tbody");
    new_output_tbody.id = "output_tbody"
    old_output_tbody = document.getElementById("output_tbody");
    old_output_tbody.parentNode.replaceChild(new_output_tbody, old_output_tbody)
    spicy_row = true;
    for (row of data) {
        table_row = new_output_tbody.insertRow();
        table_row.className = (spicy_row ? "spicy" : "normal");
        spicy_row = !spicy_row
        for (elem of row) {
            numeric_val = parseFloat(elem);
            if (!isNaN(numeric_val)) {
                elem = numeric_val.toPrecision(6)
            }
            new_cell = table_row.insertCell()
            new_cell.innerHTML = elem;
        }
    }
}

function sort_data(data) {
    data.sort(sortFunction)
    sort_order_select = document.getElementById("sort_order_select");
    if (sort_order_select.options[sort_order_select.selectedIndex].text == "Descending") {
        data = data.reverse();
    }
    return data;
}

function sortFunction(a, b) {
    i = display_filters["sortIndex"]
    aVal = parseFloat(a[i]);
    if (!isNaN(aVal)) {
        bVal = parseFloat(b[i]);
        if (aVal === bVal) {
            return 0;
        }
        else {
            return (aVal < bVal) ? -1 : 1;
        }
    }
    else if (a[i] === b[i]) {
        return 0;
    }
    else {
        return (a[i] < b[i]) ? -1 : 1;
    }
}

function search(data, search_string) {
    search_results = [];
    for (row of data) {
        if (row[0].includes(search_string.toUpperCase())) {
            search_results.push(row)
        }
    }
    data = search_results;
    return data;
}

function set_cur_page() {
    document.getElementById("page_num_input").value = display_page_num;
}

function get_num_pages() {
    num_rows = displayed_data.length
    num_pages = Math.floor(num_rows/100) + 1
    document.getElementById("total_pages").innerHTML = " of " + num_pages.toString();
    return num_pages
}

function go_to_page() {
    page_num_input = document.getElementById("page_num_input")
    page_num = parseInt(page_num_input.value, 10);
    num_pages = get_num_pages();
    if (page_num > num_pages || page_num < 1) {
        page_num = 1;
        page_num_input.value = 1;
    }
    start_index = 100 * (page_num - 1);
    end_index = start_index + 100;
    displayed_data_num_rows = displayed_data.length;
    if (displayed_data_num_rows < end_index) {
        end_index = displayed_data_num_rows;
    }
    populate_rows(start_index, end_index);
}
