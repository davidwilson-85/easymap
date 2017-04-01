/*
This .js file simply calls a PHP file via AJAX
easymap.htm --> ajax.js --> fire-wf.php --> easymap.sh --> log.log --> read_log.php ... 

*/

// Function to communicate html with php via AJAX to trigger a workflow run 
function triggerBash() {
	//alert('Button pressed');
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			document.getElementById("bash_result").innerHTML = this.responseText;
		}
	};
	
	xmlhttp.open("GET", "wf-fire.php", true);
	xmlhttp.send();
}

// Call function readLog() every few seconds
setInterval(readLog, 1000); // Second argument is the interval in milliseconds


// Function to communicate html with php via AJAX to read a log file 
function readLog() {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
		    document.getElementById("logInfo").innerHTML = this.responseText;
		}
	};
	xmlhttp.open("GET", "read-log.php", true);
	xmlhttp.send();
}



function showHint(str) {
	 if (str.length == 0) { 
	     document.getElementById("txtHint").innerHTML = "";
	     return;
	 } else {
	     var xmlhttp = new XMLHttpRequest();
	     xmlhttp.onreadystatechange = function() {
	         if (this.readyState == 4 && this.status == 200) {
	             document.getElementById("txtHint").innerHTML = this.responseText;
	         }
	     };
	     xmlhttp.open("GET", "ajax-select.php?q=" + str, true);
	     xmlhttp.send();
	 }
}
