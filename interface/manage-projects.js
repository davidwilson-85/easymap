/*
This .js file simply calls PHP files via AJAX
easymap.htm --> ajax.js --> fire-wf.php --> easymap.sh --> log.log --> read_log.php ... 

Rhis file contains the function to update project info on screen, to stop running projects,
and to delete from disk finished/killed/error projects


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

// Call function so info is shown when page is loaded
projectsInfo()

// Call function projectsInfo() every 50 seconds
setInterval(projectsInfo, 50000);


// Function to communicate html with php via AJAX to stop a project 
function stopProject(projectName) {
	if (confirm('Are you sure you want to stop running this project?')) {
		//alert('Button pressed for project ' + projectName);
		var xmlhttp = new XMLHttpRequest();
		xmlhttp.onreadystatechange = function() {
			if (this.readyState == 4 && this.status == 200) {
				// Update projects info on screen
				projectsInfo()
			}
		};
		xmlhttp.open("GET", "stop-project.php?p="+projectName, true);
		xmlhttp.send();
	}
}

// Function to communicate html with php via AJAX to remove a project from disk 
function removeProject(projectName) {
	if (confirm('Are you sure you want to permanently remove this project from disk? (Your input data will not be removed)')) {
		var xmlhttp = new XMLHttpRequest();
		xmlhttp.onreadystatechange = function() {
			if (this.readyState == 4 && this.status == 200) {
				// Update projects info on screen
				projectsInfo()
			}
		};
		xmlhttp.open("GET", "remove-project.php?p="+projectName, true);
		xmlhttp.send();
	}
}



