<?php

// This name has to be dynamic (created with javascript? / using ajax with the same php file?)
$project = 'project_name';
$loc = $project .'/log.log';

$log_contents = fopen($loc, "r");
while (!feof($log_contents)) {
	$line = fgets($log_contents);
	echo 'Line: '. $line .'<br>';
}

fclose($log_contents);

?>
