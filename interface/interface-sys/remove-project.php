<?php

$project = $_GET['project'];



$command = 'rm -rf --recursive user_projects/'. $project;

shell_exec($command);

?>
