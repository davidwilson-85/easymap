<?php

$loc = 'log.log';

$log_contents = fopen($loc, "r");
while (!feof($log_contents)) {
	$line = fgets($log_contents);
	echo 'Line: '. $line .'<br>';
}

fclose($log_contents);

?>
