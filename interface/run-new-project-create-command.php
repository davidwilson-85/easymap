<?php

/*
This file creates the command with info from javascript (user interface) and tells the shell to run
easymap.sh (the master .sh program to run easymap workflows)

Warning: If easymap.sh is run from php, the user is www-data (the apache server). Therefore www-data
must have permission to read and write in the appropriate folders.
*/

// Handle client data in JSON format
header("Content-Type: application/json");

// build a PHP variable from JSON sent using POST method
//$cmdArray = json_decode(stripslashes(file_get_contents("php://input")));
$cmdArray = json_decode(file_get_contents("php://input"));

// Elaborate the command string
//$refSeqsString = implode("+", $cmdArray[4]);
//$cmdArray[4] = $refSeqsString;
$cmdString = implode(" ", $cmdArray);

// Run workflow
shell_exec('cd ..; '. $cmdString);

// Only for development
//$file = fopen('file.txt', 'w');
//fwrite($file, $cmdString);
//fclose($file);

//echo 'success';

?>
