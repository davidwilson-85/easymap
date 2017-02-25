
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
	
function verifyProjectName(){
	var text = document.getElementById("form1").projectName.value;
	if( /[^a-zA-Z0-9]/.test( text ) ) {
		alert('Input is not alphanumeric');
	} else if (text == '') {
		alert('Yoy must give a name to the project');
	} else {
		cmdArgs[0] = document.getElementById("form1").projectName.value;
		updateCmd()
	}
}

	function projectName() {
		cmdArgs[0] = document.getElementById("form1").projectName.value;
		updateCmd()
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
			cmdArgs[1] = 'ins';
			document.getElementById("expDataIns").style.display = "inline";
			document.getElementById("expDataSnp").style.display = "none";
		} else {
			cmdArgs[1] = 'snp';
			document.getElementById("expDataIns").style.display = "none";
			document.getElementById("expDataSnp").style.display = "inline";
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
			cmdArgs[2] = 'exp';
			document.getElementById("expDataInterface").style.display = "inline";
			document.getElementById("simDataInterface").style.display = "none";
		} else {
			cmdArgs[2] = 'sim';
			document.getElementById("expDataInterface").style.display = "none";
			document.getElementById("simDataInterface").style.display = "inline";
		}
		updateCmd();
	}
	
	// Determine option button selected and define the appropriate command argument
	function buttons_libType() {
		var options = document.getElementsByClassName("libType");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button5') {
			cmdArgs[3] = 'se';
			document.getElementById("expDataSingle").style.display = "inline";
			document.getElementById("expDataPaired").style.display = "none";
			document.getElementById("expDataSingleTwosamples").style.display = "inline";
			document.getElementById("expDataPairedTwosamples").style.display = "none";
		} else {
			cmdArgs[3] = 'pe';
			document.getElementById("expDataSingle").style.display = "none";
			document.getElementById("expDataPaired").style.display = "inline";
			document.getElementById("expDataSingleTwosamples").style.display = "none";
			document.getElementById("expDataPairedTwosamples").style.display = "inline";
		}
		updateCmd();
	}
	
	// Determine all the file names selected, add them to array, and then to argument
	function refSeqs() {
		var contigs = document.getElementById("refSeqs");
		var contigsList = [];
		
		for (var i=0; i<contigs.length; i++) {
			if (contigs[i].selected == true) {
				contigsList.push(contigs[i].value);
			}
		}
		cmdArgs[4] = contigsList;
		updateCmd();
	}
	
	function processSingleSelectors() {		
		if (this.id == 'insSeq') {
			cmdArgs[5] = this.value;
		}
		if (this.id == 'gffFile') {
			cmdArgs[9] = this.value;
		}
		if (this.id == 'annFile') {
			cmdArgs[10] = this.value;
		}
		if (this.id == 'readsSingle') {
			cmdArgs[6] = this.value;
		}
		if (this.id == 'readsForward') {
			cmdArgs[7] = this.value;
		}
		if (this.id == 'readsReverse') {
			cmdArgs[8] = this.value;
		}
		
		updateCmd();	
	}
	
	// End of functions
	
	
	// Define array with all the command arguments
	var cmdArgs = ['project', 'workflow', 'dataSource', 'libType', 'refSeqs', 'insSeq',
						'readsS', 'readsF', 'readsR', 'gffFile', 'annFile', 'simMut', 'simRecsel',
						'simSeq'];
	
	// React to interactions with project name field
	document.getElementById("form1").projectName.onfocus = resetTextField;
	document.getElementById("form1").projectName.onblur = verifyProjectName;
	//document.getElementById("form1").projectName.onblur = projectName;
	
	// Change the content of the attribute href of a specific link
	updateCmd();
	
	// React to interactions with 2-way selectors
	document.getElementById("button1").onclick = buttons_analysisType;
	document.getElementById("button2").onclick = buttons_analysisType;
	document.getElementById("button3").onclick = buttons_dataSource;
	document.getElementById("button4").onclick = buttons_dataSource;
	document.getElementById("button5").onclick = buttons_libType;
	document.getElementById("button6").onclick = buttons_libType;
	
	// React to interactions with genome contigs selector
	document.getElementById("form1").refSeqs.onblur = refSeqs;
	
	// React to single selectots (insSeq, gffFile, annFile...)
	document.getElementById("form1").insSeq.onblur = processSingleSelectors;
	document.getElementById("form1").gffFile.onblur = processSingleSelectors;
	document.getElementById("form1").annFile.onblur = processSingleSelectors;
	
	document.getElementById("form1").readsSingle.onblur = processSingleSelectors;
	document.getElementById("form1").readsForward.onblur = processSingleSelectors;
	document.getElementById("form1").readsReverse.onblur = processSingleSelectors;

}






