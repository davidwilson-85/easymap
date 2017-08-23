/*
This .js file simply calls PHP files via AJAX
easymap.htm --> ajax.js --> fire-wf.php --> easymap.sh --> log.log --> read_log.php ... 

*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// THIS SECTION DEALS WITH GENERAL FUNCTIONALITY OF THE PAGE, SUSH AS CHECKING CONFIG FILE OR TRIGGERING A PROJECT
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to communicate html with php via AJAX to check if user_projects directory
// is over config/config>user_projects-size-limit
// This function has to be triggered whe page loads
function allowNewProject() {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			var response = this.responseText;
			var responseFields = response.split(",");
			if (responseFields[0] >= 100) {
				document.getElementById("runNewProject").style.display = "none";
				document.getElementById("sizeWarning").style.display = "block";
				document.getElementById("sizeWarning").innerHTML = "WARNING: You are over the maximum space limit set by the machine administrator (" + responseFields[2] + " Gb). Delete past projects to free disk space.";
			} else if (responseFields[0] >= 75 && responseFields[0] < 100) {
				document.getElementById("sizeWarning").style.display = "block";
				document.getElementById("sizeWarning").innerHTML = "WARNING: You are over 75% of the maximum space limit set by the machine administrator (" + responseFields[2] + " Gb).";
			}
			if (responseFields[1] >= 100) {
				document.getElementById("runNewProject").style.display = "none";
				document.getElementById("simultWarning").style.display = "block";
				document.getElementById("simultWarning").innerHTML = "WARNING: You have reached the maximum number of simultaneously running projects set by the machine administrator (" + responseFields[3] + " projects). Wait until currently running jobs finish.";
			}
		}
	};
	xmlhttp.open("GET", "allow-new-project.php", true);
	xmlhttp.send();
}

// Function to list all genome reference files to be displayed inside <select multiple> tag
// This function has to be triggered whe page loads
function listInputFiles() {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {

			// Parse response in JSON format
			var inputFilesresponse = JSON.parse(this.responseText);

			// Create interface to select reference sequence contigs and insertion sequence
			var fastaFiles = inputFilesresponse[0];
			var refFiles = document.getElementById('refFileSelector');
			var insFiles = document.getElementById('insFileSelector');
			insFiles.options[insFiles.options.length] = new Option('Select a file', 'XXX');
			for (i = 0; i < fastaFiles.length; i++) {
				refFiles.options[refFiles.options.length] = new Option(fastaFiles[i], fastaFiles[i]);
				insFiles.options[insFiles.options.length] = new Option(fastaFiles[i], fastaFiles[i]);
			}

			// Create interface to select GFF and ANN files
			var otherFiles = inputFilesresponse[2];
			var gffFiles = document.getElementById('gffFileSelector');
			var annFiles = document.getElementById('annFileSelector');
			gffFiles.options[gffFiles.options.length] = new Option('Select a file', 'XXX');
			annFiles.options[annFiles.options.length] = new Option('Select a file', 'XXX');
			for (i = 0; i < otherFiles.length; i++) {
				gffFiles.options[gffFiles.options.length] = new Option(otherFiles[i], otherFiles[i]);
				annFiles.options[annFiles.options.length] = new Option(otherFiles[i], otherFiles[i]);
			}

			// Create interfaces to select fastq files
			var fastqFiles = inputFilesresponse[1];
			var readsProblem = document.getElementById('readsProblemSelector');
			var readsControl = document.getElementById('readsControlSelector');
			for (i = 0; i < fastqFiles.length; i++) {
				readsProblem.options[readsProblem.options.length] = new Option(fastqFiles[i], fastqFiles[i]);
				readsControl.options[readsControl.options.length] = new Option(fastqFiles[i], fastqFiles[i]);
			}
		}
	};
	xmlhttp.open("GET", "run-new-project-listInputFiles.php", true);
	xmlhttp.send();
}

// Function to communicate html with php via AJAX to trigger a workflow run 
function runProject() {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			document.getElementById("commandResponse").innerHTML = this.responseText;
		}
	};
	var arguments = "test-arg";
	xmlhttp.open("GET", "run-new-project-create-command.php?args=" + arguments, true);
	xmlhttp.send();
}

// Trigger required functions when page loads
allowNewProject();
listInputFiles()


///////////////////////////////////////////////////////////////////////////////////////////////////////
// THIS SECTION DEALS WITH THE DYNAMIC FORMATTING OF THE FORM AND WITH FIELDS VALIDATION
///////////////////////////////////////////////////////////////////////////////////////////////////////


function resetTextField() {
	this.value = "";
}
	
window.onload = function() {
	
	// Functions that need to be declared after document has completely loaded
	
	// Uptdate the command arguments in the screen
	// This function is called in many other functions after updating the value
	// of an arguments.
	function updateCmd() {
			document.getElementById("commandString").innerHTML = cmdArgs;
	}
	
	//
	function verifyProjectName(){
		var text = document.getElementById("form1").projectName.value;
		if(/[^a-zA-Z0-9]/.test( text) ) {
			//alert('Input is not alphanumeric');
			projectNameValidationInfoMessage = 'Input is not alphanumeric';
			document.getElementById("projectNameValidationInfo").innerHTML = projectNameValidationInfoMessage;
			document.getElementById("projectNameValidationInfo").style.display = "block";
		} else if (text == '') {
			//alert('Yoy must give a name to the project');
			projectNameValidationInfoMessage = 'You must give a name to the project';
			document.getElementById("projectNameValidationInfo").innerHTML = projectNameValidationInfoMessage;
			document.getElementById("projectNameValidationInfo").style.display = "block";
		} else {
			cmdArgs[1] = document.getElementById("form1").projectName.value;
			updateCmd()
			projectNameValidationInfoMessage = '';
			document.getElementById("projectNameValidationInfo").innerHTML = projectNameValidationInfoMessage;
			document.getElementById("projectNameValidationInfo").style.display = "none";
		}
	}
	
	// Determine option button selected and define the appropriate command argument
	function buttons_analysisType() {
		var options = document.getElementsByClassName("analysisType");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button1') {
			cmdArgs[2] = 'ins';
			document.getElementById("insSeqField").style.display = "inline";			//////////////////////////
			document.getElementById("readsControl").style.display = "none";				//////////////////////////
			document.getElementById("backgroundCrossCtype").style.display = "none";
			//document.getElementById("simDataIns").style.display = "inline";			//////////////////////////
			//document.getElementById("simDataSnp").style.display = "none";				//////////////////////////
		} else {
			cmdArgs[2] = 'snp';
			document.getElementById("insSeqField").style.display = "none";
			document.getElementById("readsControl").style.display = "inline";
			document.getElementById("backgroundCrossCtype").style.display = "block";
			//document.getElementById("simDataIns").style.display = "none";
			//document.getElementById("simDataSnp").style.display = "inline";
		}
		updateCmd();
	}
	
	// Determine option button selected and define the appropriate command argument
	function buttons_dataSource() {
		var options = document.getElementsByClassName("dataSource");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button3') {
			cmdArgs[3] = 'exp';
			document.getElementById("expDataInterface").style.display = "inline";
			document.getElementById("simDataInterface").style.display = "none";
		} else {
			cmdArgs[3] = 'sim';
			document.getElementById("expDataInterface").style.display = "none";
			document.getElementById("simDataInterface").style.display = "inline";
		}
		updateCmd();
	}
	
	// Determine all the reference file names selected, add them to array, and then to command argument
	function checkRefSeqs() {
		var contigs = document.getElementById("refFileSelector");
		var contigsList = [];
		
		for (var i=0; i<contigs.length; i++) {
			if (contigs[i].selected == true) {
				contigsList.push(contigs[i].value);
			}
		}
		cmdArgs[4] = contigsList;
		updateCmd();
	}

	// Update command arguments after each user interaction with sinlge selectors
	function processSingleSelectors() {
		if (this.id == 'insFileSelector') {
			cmdArgs[5] = this.value;
		}
		if (this.id == 'gffFileSelector') {
			cmdArgs[6] = this.value;
		}
		if (this.id == 'annFileSelector') {
			cmdArgs[7] = this.value;
		}
		updateCmd();	
	}

	// Mutant background: determine option button selected and define the appropriate command argument
	function buttons_mutBackground() {
		var options = document.getElementsByClassName("mutBackground");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button11') {
			cmdArgs[16] = 'ref';
		} else {
			cmdArgs[16] = 'noref';
		}
		updateCmd();
	}

	// Mapping cross preformed: determine option button selected and define the appropriate command argument
	function buttons_crossType() {
		var options = document.getElementsByClassName("crossType");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button13') {
			cmdArgs[17] = 'bc';
		} else {
			cmdArgs[17] = 'oc';
		}
		updateCmd();
	}

	// Origin of the control reads: determine option button selected and define the appropriate command argument
	function buttons_contType() {
		var options = document.getElementsByClassName("contType");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button15') {
			cmdArgs[18] = 'par';
			cmdArgs[19] = 'mut';
		} else if (checkedOption == 'button16') {
			cmdArgs[18] = 'par';
			cmdArgs[19] = 'nomut';
		} else {
			cmdArgs[18] = 'f2wt';
			cmdArgs[19] = 'n/p';
		}
		updateCmd();
	}

	// Check if combination of mutat background, cross performed, and origin of control reads, is supported
	function checkBackgroundCrossCtypeIntermediateCheck() {
		if (cmdArgs[16] == "ref" && cmdArgs[17] == "oc") {
			document.getElementById("backgroundCrossCtypeWarnMsg").style.display = "block";
		} else {
			document.getElementById("backgroundCrossCtypeWarnMsg").style.display = "none";
		}
	}

/*
	// Check if combination of mutant background, cross performed, and origin of control reads, is supported
	function checkBackgroundCrossCtypeFinalCheck() {
		if (cmdArgs[16] == "n/p" || cmdArgs[17] == "n/p" || cmdArgs[18] == "n/p") {
			document.getElementById("backgroundCrossCtypeWarnMsg").style.display = "block";
		} else {
			document.getElementById("backgroundCrossCtypeWarnMsg").style.display = "none";
		}
	}
*/
	// Check reads selectors (max two files selected per sample) and update command argument(s)
	function checkProblemReads() {
		var reads = document.getElementById("readsProblemSelector");
		var readsList = [];
		
		for (var i=0; i<reads.length; i++) {
			if (reads[i].selected == true) {
				readsList.push(reads[i].value);
			}
		}

		if (readsList.length == 1) {
			cmdArgs[8] = readsList; cmdArgs[9] = 'XXX'; cmdArgs[10] = 'XXX'; cmdArgs[11] = 'se';
			updateCmd();
			// Hide error message
			document.getElementById("readsProblemWarnMsg").style.display = "none";
		} else if (readsList.length == 2) {
			cmdArgs[8] = 'XXX'; cmdArgs[9] = readsList[0]; cmdArgs[10] = readsList[1]; cmdArgs[11] = 'pe';
			updateCmd();
			// Hide error message
			document.getElementById("readsProblemWarnMsg").style.display = "none";
		} else {
			cmdArgs[8] = 'XXX'; cmdArgs[9] = 'XXX'; cmdArgs[10] = 'XXX'; cmdArgs[11] = 'XXX';
			updateCmd();
			// Display error message
			document.getElementById("readsProblemWarnMsg").style.display = "block";
		}
	}

	function checkControlReads() {
		var reads = document.getElementById("readsControlSelector");
		var readsList = [];
		
		for (var i=0; i<reads.length; i++) {
			if (reads[i].selected == true) {
				readsList.push(reads[i].value);
			}
		}

		if (readsList.length == 1) {
			cmdArgs[12] = readsList; cmdArgs[13] = 'XXX'; cmdArgs[14] = 'XXX'; cmdArgs[15] = 'se';
			updateCmd();
			// Hide error message
			document.getElementById("readsControlWarnMsg").style.display = "none";
		} else if (readsList.length == 2) {
			cmdArgs[12] = 'XXX'; cmdArgs[13] = readsList[0]; cmdArgs[14] = readsList[1]; cmdArgs[15] = 'pe';
			updateCmd();
			// Hide error message
			document.getElementById("readsControlWarnMsg").style.display = "none";
		} else {
			cmdArgs[12] = 'XXX'; cmdArgs[13] = 'XXX'; cmdArgs[14] = 'XXX'; cmdArgs[15] = 'XXX';
			updateCmd();
			// Show error message
			document.getElementById("readsControlWarnMsg").style.display = "block";
		}
	}
	
	// End of functions
	
	
	// Define array with all the command arguments
	var cmdArgs = ['./easymap.sh','project_name','workflow','data_source','ref_seq','ins_seq','gff_file',
					'ann_file','read_s','read_f','read_r','lib_type_sample','read_s_ctrl','read_f_ctrl',
					'read_r_ctrl','lib_type_ctrl','n/p','n/p','n/p','n/p','sim_mut','sim_recsel','sim_seq'];

	// Create the command string for the first time (for development purposes only)
	updateCmd();
	
	// React to interactions with text inputs
	// Reset default content when user clicks on input bux
	document.getElementById("form1").projectName.onfocus = resetTextField;
	document.getElementById("form1").simMutNbr.onfocus = resetTextField;
	document.getElementById("form1").simSeqRD.onfocus = resetTextField;
	document.getElementById("form1").simSeqRL.onfocus = resetTextField;
	document.getElementById("form1").simSeqRD.onfocus = resetTextField;
	document.getElementById("form1").simSeqBER.onfocus = resetTextField;
	document.getElementById("form1").simSeqGBS.onfocus = resetTextField;
	
	// Verify input of text fields
	document.getElementById("form1").projectName.onblur = verifyProjectName;
	
	// React to interactions with the main 2-way selectors
	document.getElementById("button1").onclick = buttons_analysisType;
	document.getElementById("button2").onclick = buttons_analysisType;
	document.getElementById("button3").onclick = buttons_dataSource;
	document.getElementById("button4").onclick = buttons_dataSource;
	
	// React to interactions with genome contigs selector
	document.getElementById("form1").refFileSelector.onblur = checkRefSeqs;
	
	// React to single selectors (insSeq, gffFile, annFile...)
	document.getElementById("form1").insFileSelector.onblur = processSingleSelectors;
	document.getElementById("form1").gffFileSelector.onblur = processSingleSelectors;
	document.getElementById("form1").annFileSelector.onblur = processSingleSelectors;
	
	// React to interactions with the MbS 2-way selectors
	document.getElementById("button11").onclick = buttons_mutBackground;
	document.getElementById("button12").onclick = buttons_mutBackground;
	document.getElementById("button13").onclick = buttons_crossType;
	document.getElementById("button14").onclick = buttons_crossType;
	document.getElementById("button15").onclick = buttons_contType;
	document.getElementById("button16").onclick = buttons_contType;
	document.getElementById("button17").onclick = buttons_contType;

	document.getElementById("button17").onmouseover = checkBackgroundCrossCtypeIntermediateCheck;


	//React to interactions with reads selectors
	document.getElementById("form1").readsProblemSelector.onblur = checkProblemReads;
	document.getElementById("form1").readsControlSelector.onblur = checkControlReads;

	// React to interactions with backgroundCrossCtype buttons
	document.getElementById("backgroundCrossCtype").onmouseout = checkBackgroundCrossCtypeIntermediateCheck;
}

// TO DO:
// Add verification to check that all the three 2-way buttons have been clicked
// When user clicks on submit, check that all validations passed. This could be done with
// a flag variable.
//
//
//
//
//
//
//
//
//
//





