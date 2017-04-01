<?php

$project = $_POST['project'];

$project = '';



$command = 'rm -rf --recursive '. $project;

$output = shell_exec($command);

?>
