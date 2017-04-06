/*
This .js file simply calls PHP files via AJAX
easymap.htm --> ajax.js --> fire-wf.php --> easymap.sh --> log.log --> read_log.php ... 

*/

// Enclose all code inside a onload condition?????????????????????????


// Function to communicate html with php via AJAX to retrieve projects info 
function projectsInfo() {
	//alert('Button pressed');
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			document.getElementById("projectsInfo").innerHTML = this.responseText;
		}
	};
	
	xmlhttp.open("GET", "manage-projects.php", true);
	xmlhttp.send();
}

// Call function so
projectsInfo()

// Call function projectsInfo() every 50 seconds
setInterval(projectsInfo, 50000);

