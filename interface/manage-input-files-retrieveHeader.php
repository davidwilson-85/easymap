<?php

$file = $_GET['f'];

$nbrOfLines = 100;

$command = 'head -'. $nbrOfLines .' ../user_data/'. $file;

$preview = shell_exec($command);

echo '<p>File: '. $file .'.</p><p>(showing first '. $nbrOfLines .' lines)</p>';

echo '<pre>'. $preview .'</pre>';


?>
