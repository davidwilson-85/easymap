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


// Define a string that contais all the arguments
$arguments = $project_name .' '. $workflow .' '. $data_source;

// Run workflow
$command = './test-workflow.sh '. $arguments;
$output = shell_exec($command);



?>
