<?php
// Start of code for dealing with password

session_start();
$reference_password = $_SESSION['reference_password'];

if (isset($_POST['password'])) {
	$_SESSION['password'] = $_POST['password'];
	$password = $_SESSION['password'];
} else {
	@$password = $_SESSION['password'];
}

if ($password != $reference_password) {

	echo '
		<!DOCTYPE html>
		<html>
			<head>
				<title></title>
			</head>
			<body>
				Wrong password. <a href="index.php">Try again</a>
			</body>
		</html>
	';

}

if ($password == $reference_password) {
// End of code for dealing with password
?>

<!DOCTYPE html>
<html>
<title>Easymap</title>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">

<link rel="stylesheet" href="w3c.css">
<link rel="stylesheet" href="style.css">
<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Poppins">
<link rel="stylesheet" type="text/css" href="run-new-project-design-project.css">

<script type="text/javascript" src="run-new-project.js"></script>

<style>
body,h1,h2,h3,h4,h5 {font-family: "Poppins", sans-serif}
body {font-size:16px;}
.w3-half img{margin-bottom:-6px;margin-top:16px;opacity:0.8;cursor:pointer}
.w3-half img:hover{opacity:1}
</style>

<body>

<!-- Sidebar/menu -->
<nav class="w3-sidebar w3-red w3-collapse w3-top w3-large w3-padding" style="z-index:3;width:300px;font-weight:bold;" id="mySidebar"><br>
  <a href="javascript:void(0)" onclick="w3_close()" class="w3-button w3-hide-large w3-display-topleft" style="width:100%;font-size:22px">Close Menu</a>
  <div class="w3-container">
    <h3 class="w3-padding-64"><b>Easymap</b></h3>
  </div>
  <div class="w3-bar-block">
    <a href="manage-input-files.php" onclick="w3_close()" class="w3-bar-item w3-button w3-hover-white">Manage input files</a> 
    <a href="manage-projects.php" onclick="w3_close()" class="w3-bar-item w3-button w3-hover-white">Manage projects</a> 
    <a href="run-new-project.php" onclick="w3_close()" class="w3-bar-item w3-button w3-hover-white">Run new project</a> 
    <a href="documentation.php" onclick="w3_close()" class="w3-bar-item w3-button w3-hover-white">Documentation</a>
  </div>
</nav>

<!-- Top menu on small screens -->
<header class="w3-container w3-top w3-hide-large w3-red w3-xlarge w3-padding">
  <a href="javascript:void(0)" class="w3-button w3-red w3-margin-right" onclick="w3_open()">☰</a>
  <span>Easymap</span>
</header>

<!-- Overlay effect when opening sidebar on small screens -->
<div class="w3-overlay w3-hide-large" onclick="w3_close()" style="cursor:pointer" title="close side menu" id="myOverlay"></div>

<!-- beginning of page content -->
<div class="w3-main" style="margin-left:340px;margin-right:40px">

  <div class="w3-container" style="margin-top:75px">
    <h1 class="w3-xxxlarge w3-text-red"><b>Run new project</b></h1>
    <hr style="width:50px;border:5px solid red" class="w3-round">
    <h3>Design and execute a new project</h3>
    
    <!-- If the size of data or the number or currently running projects is over the limits,
    display Warning messages and do not display id runNewProject -->
    
    <div id="sizeWarning" class="warningMessage"></div>
    <div id="simultWarning" class="warningMessage"></div>
    
    <div id="runNewProject">
			
		<form id="form1">
			
			<hr style="width: 100%; border: 2px solid rgb(150,150,150)" class="w3-round">
			
			<p>
				Give a name to this project:<br>
				<input type="text" name="projectName" value="My project" />
				<div id="projectNameValidationInfo" class="warningMessage"></div>
			</p>
			
			Mapping-by-sequencing strategy:
			<div class="buttons-wrap">
				<div class="mx-button">
					<input type="radio" class="analysisType" name="mx12" id="button1" />
					<label for="button1" unselectable>Tagged sequence mapping</label>
				</div>
				<div class="mx-button">
					<input type="radio" class="analysisType" name="mx12" id="button2" />
					<label for="button2" unselectable>Linkage analysis mapping</label>
				</div>
				<div class="clear-floats"></div>
			</div>
			
			Data source:
			<div class="buttons-wrap">
				<div class="mx-button">
					<input type="radio" class="dataSource" name="mx34" id="button3" />
					<label for="button3" unselectable>Use my own data</label>
				</div>
				<div class="mx-button">
					<input type="radio" class="dataSource" name="mx34" id="button4" />
					<label for="button4" unselectable>Simulate data</label>
				</div>
				<div class="clear-floats"></div>
			</div>
	
			Type of NGS library:
			<div class="buttons-wrap">
				<div class="mx-button">
					<input type="radio" class="libType" name="mx56" id="button5" />
					<label for="button5" unselectable>Single end</label>
				</div>
				<div class="mx-button">
					<input type="radio" class="libType" name="mx56" id="button6" />
					<label for="button6" unselectable>Paired end</label>
				</div>
				<div class="clear-floats"></div>
			</div>
			
			<hr style="width: 100%; border: 2px solid rgb(150,150,150)" class="w3-round">
			
			Reference sequence (You can select multiple files by pressing and holding the Ctrl/Cmd key):<br>
			<div id="refSeqs"></div>
			
			<div id="insSeqField"> <!-- Not displayed by default -->
				Insertion sequence file:
				<div id="insSeq"></div>
			</div>
			
			GFF3 file:<br>
			<div id="gffFile"></div>

			Gene functional annotation file:<br>
			<div id="annFile"></div>

		
			<div id="expDataInterface">
				Experimental data chosen
						
				<div id="expDataIns">
					<div id="expDataSingle">
						<select id="readsSingle">
							<option value="default">Choose single-end reads</option>
							<option value="insertion.fa">insertion.fa</option>
							<option value="arabidopsis_gff.gff">arabidopsis_gff.gff</option>
							<option value="at-gene-funct-tair10.txt">at-gene-funct-tair10.txt</option>
							<option value="reads_f.fq">reads_f.fq</option>
							<option value="reads_r.fq">reads_r.fq</option>
							<option value="reads_se.fq">reads_se.fq</option>
						</select>
					</div>
				
					<div id="expDataPaired">
						<select id="readsForward">
							<option value="default">Choose forward reads</option>
							<option value="insertion.fa">insertion.fa</option>
							<option value="arabidopsis_gff.gff">arabidopsis_gff.gff</option>
							<option value="at-gene-funct-tair10.txt">at-gene-funct-tair10.txt</option>
							<option value="reads_f.fq">reads_f.fq</option>
							<option value="reads_r.fq">reads_r.fq</option>
							<option value="reads_se.fq">reads_se.fq</option>
						</select>
						<select id="readsReverse">
							<option value="default">Choose reverse reads</option>
							<option value="insertion.fa">insertion.fa</option>
							<option value="arabidopsis_gff.gff">arabidopsis_gff.gff</option>
							<option value="at-gene-funct-tair10.txt">at-gene-funct-tair10.txt</option>
							<option value="reads_f.fq">reads_f.fq</option>
							<option value="reads_r.fq">reads_r.fq</option>
							<option value="reads_se.fq">reads_se.fq</option>
						</select>
					</div>
				</div>
				
				<div id="expDataSnp">
				
					<div id="expDataSnpInfo">				
						<div>
							Mutant strain<br>
							<input type="radio" id="refStrainYes" name="isRefStrain">
							<label for="refStrainYes">The same as reference strain</label>
							<br>
							<input type="radio" id="refStrainNo" name="isRefStrain">
							<label for="refStrainNo">Other than the reference strain</label>
						</div>
						<br>
						<div>
							Type of mapping cross performed<br>
							<input type="radio" id="crossTypeBack" name="crossType">
							<label for="crossTypeBack">Backcross</label>
							<br>
							<input type="radio" id="crossTypeOut" name="crossType">
							<label for="crossTypeOut">Outcross</label>
						</div>
						<br>
						<div>
							Parental reads provided<br>
							<input type="radio" id="parentalReadsMutant" name="parentalReads">
							<label for="parentalReadsMutant">The mutant's</label>
							<br>
							<input type="radio" id="parentalReadsNoMutant" name="parentalReads">
							<label for="parentalReadsNoMutant">The other</label>
						</div>
					</div>
				
					<div id="expDataSingleTwosamples">
						<select id="readsSinglePar">
							<option value="default">Parental - choose reads</option>
							<option value="insertion.fa">insertion.fa</option>
							<option value="arabidopsis_gff.gff">arabidopsis_gff.gff</option>
							<option value="at-gene-funct-tair10.txt">at-gene-funct-tair10.txt</option>
							<option value="reads_f.fq">reads_f.fq</option>
							<option value="reads_r.fq">reads_r.fq</option>
							<option value="reads_se.fq">reads_se.fq</option>
						</select>
						<br>
						<select id="readsSingleF2">
							<option value="default">F2 - choose reads</option>
							<option value="insertion.fa">insertion.fa</option>
							<option value="arabidopsis_gff.gff">arabidopsis_gff.gff</option>
							<option value="at-gene-funct-tair10.txt">at-gene-funct-tair10.txt</option>
							<option value="reads_f.fq">reads_f.fq</option>
							<option value="reads_r.fq">reads_r.fq</option>
							<option value="reads_se.fq">reads_se.fq</option>
						</select>
					</div>
			
					<div id="expDataPairedTwosamples">
						<select id="readsForwardPar">
							<option value="default">Parental - choose forward reads</option>
							<option value="insertion.fa">insertion.fa</option>
							<option value="arabidopsis_gff.gff">arabidopsis_gff.gff</option>
							<option value="at-gene-funct-tair10.txt">at-gene-funct-tair10.txt</option>
							<option value="reads_f.fq">reads_f.fq</option>
							<option value="reads_r.fq">reads_r.fq</option>
							<option value="reads_se.fq">reads_se.fq</option>
						</select>
						<select id="readsReversePar">
							<option value="default">Parental - choose reverse reads</option>
							<option value="insertion.fa">insertion.fa</option>
							<option value="arabidopsis_gff.gff">arabidopsis_gff.gff</option>
							<option value="at-gene-funct-tair10.txt">at-gene-funct-tair10.txt</option>
							<option value="reads_f.fq">reads_f.fq</option>
							<option value="reads_r.fq">reads_r.fq</option>
							<option value="reads_se.fq">reads_se.fq</option>
						</select>
						<br>
						<select id="readsReversePar">
							<option value="default">F2 - choose forward reads</option>
							<option value="insertion.fa">insertion.fa</option>
							<option value="arabidopsis_gff.gff">arabidopsis_gff.gff</option>
							<option value="at-gene-funct-tair10.txt">at-gene-funct-tair10.txt</option>
							<option value="reads_f.fq">reads_f.fq</option>
							<option value="reads_r.fq">reads_r.fq</option>
							<option value="reads_se.fq">reads_se.fq</option>
						</select>
						<select id="readsReverseF2">
							<option value="default">F2 - choose reverse reads</option>
							<option value="insertion.fa">insertion.fa</option>
							<option value="arabidopsis_gff.gff">arabidopsis_gff.gff</option>
							<option value="at-gene-funct-tair10.txt">at-gene-funct-tair10.txt</option>
							<option value="reads_f.fq">reads_f.fq</option>
							<option value="reads_r.fq">reads_r.fq</option>
							<option value="reads_se.fq">reads_se.fq</option>
						</select>
					</div>
				</div>
			</div>
		
			<div id="simDataInterface">
				Simulated data chosen
				<div id="simDataIns">				
					<div id="simMutInterface">
						Number of mutations:<br>
						<input type="text" name="simMutNbr" class="simNumericInput" value="" />
						<br><br>
						Mode:<br>
						<input type="radio" name="simMutMode" id="simMutModeD">
						<label for="simMutModeD">Drift</label>
						<br>
						<input type="radio" name="simMutMode" id="simMutModeE">
						<label for="simMutModeE">EMS (GC > AT)</label>
						<br>
						<input type="radio" name="simMutMode" id="simMutModeL">
						<label for="simMutModeL">Large DNA insertions</label>
					</div>
			
					<div id="simRecselInterface">
						<br><br>No RecSel here<br><br>
					</div>
			
					<div id="simSeqInterface">
						Read depth (X):<br>
						<input type="text" name="simSeqRD" class="simNumericInput" value="40" />
						<br>
						Read length [mean, SD] (nt):<br>
						<input type="text" name="simSeqRL" class="simNumericInput" value="90,0" />
						<br>
						<div id="simSeqFL">
							Fragment length [mean, SD] (nt):<br>
							<input type="text" name="simSeqFL" class="simNumericInput" value="500,100" />
						</div>
						Basecalling error rate (%):<br>
						<input type="text" name="simSeqBER" class="simNumericInput" value="1" />
						<br>
						GC bias strength (%):<br>
						<input type="text" name="simSeqGBS" class="simNumericInput" value="50" />
						<br>
					</div>
			
				</div>
				<div id="simDataSnp">
			
					Simulate Mapping by Sequencing data
			
				</div>
			
			</div>
	
			<div id="formButtons">
				<input type="button" class="button" id="runButton" value="Run analysis" /> 
				<input type="button" class="button" id="otherButton" value="Other button" />
			</div>
		
		</form>
		
		<br><br><br>

		<div id="command">
			<p>
				$project_name $workflow $data_source $lib_type $ref_seq $ins_seq $read_s $reads_f $reads_r $gff_file $ann_file $sim-mut $sim-recsel $sim-seq
			</p>
			<p id="commandString"></p>
		</div>	
    	
    	<a href="manage-projects.php" class="button" onclick="runProject()">Run workflow</a>
    </div>
    
  </div>
  

<!-- End of page content -->
</div>

<!-- W3.CSS Container -->
<div class="w3-light-grey w3-container w3-padding-32" style="margin-top:75px;padding-right:58px"><p class="w3-right">Template: w3c</a></p></div>

<script>
// Script to open and close sidebar
function w3_open() {
    document.getElementById("mySidebar").style.display = "block";
    document.getElementById("myOverlay").style.display = "block";
}
 
function w3_close() {
    document.getElementById("mySidebar").style.display = "none";
    document.getElementById("myOverlay").style.display = "none";
}
/*
// Modal Image Gallery
function onClick(element) {
  document.getElementById("img01").src = element.src;
  document.getElementById("modal01").style.display = "block";
  var captionText = document.getElementById("caption");
  captionText.innerHTML = element.alt;
}
*/
</script>

</body>
</html>

<?php
// Code for dealing with password
}
?>
