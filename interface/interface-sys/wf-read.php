<?php


// This program has to be executed every X seconds via AJAX in order to keep the user screen updated


// This name has to be dynamic (created with javascript? / using ajax with the same php file?)
$project = 'project_name';
$loc = 'log.log';

$log_contents = fopen($loc, "r");
while (!feof($log_contents)) {
	$line = fgets($log_contents);
	echo 'Line: '. $line .'<br>';
}

fclose($log_contents);

?>
