<?php

/*
This file creates the command with info from javascript (user interface) and tells the shell to run
easymap.sh (the master .sh program to run easymap workflows)


Warning: If easymap.sh is run from fire_wf.php, the user is www-data (the apache server). Therefore www-data
must have permission to read and write in the appropriate folders.



*/

// Define the names of the input files provided by the user
$project_name = 'project';
$workflow = 'ins';
$data_source = 'sim';
$lib_type = 'se';


// Define a string that contais all the arguments
$arguments = $project_name .' '. $workflow .' '. $data_source .' '. $lib_type;

// Run workflow
$command = './program.sh '. $arguments;

//echo $command .'<br>';

$output = shell_exec($command);

/*
// Read workflow log file
$loc = 'log.log';
$log_contents = fopen($loc, "r");
while (!feof($log_contents)) {
	$line = fgets($log_contents);
	echo 'Line: '. $line .'<br>';
}

fclose($log_contents);


echo 'Done!';
*/
?>
