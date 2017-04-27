<?php

$location = 'log.log';

$log_contents = fopen($location, "r");
while (!feof($log_contents)) {
	$line = fgets($log_contents);
	echo 'Line: '. $line .'<br>';
}

fclose($log_contents);

?>
