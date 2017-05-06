<?php

$project_name = $_GET['p'];

$location = '../user_projects/'. $project_name .'/3_workflow_output/report.html';

// Read HTML file with report as a text file and line by line
// Send each line to HTML via AJAX
$report_contents = fopen($location, "r");
while (!feof($report_contents)) {
	$line = fgets($report_contents);
	echo $line .'<br>';
}
fclose($report_contents);

?>
