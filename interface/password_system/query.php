<?php //Code for dealing with password======================================================================

if (isset($_SESSION['password'])) {

	session_start(); 
	$password = $_SESSION['password'];
	echo 'Password is: '. $password;

} else {

	$password = $_POST['password'];
	
	session_start();
	$_SESSION['password'] = $password;

}




?>

<?php

if ($password == "elche") { ?>

<html>
<head>
<title>Correct password</title>
</head>
<body>
<p>UMH leaf mutant DB - Correct password</p>
<p><a href="page2.php">Go to page2.php</a></p>
</body>
</html>

<?php } else { ?>

<html>
<head>
<title>Wrong password</title>
</head>
<body>
<p>UMH leaf mutant DB - Wrong password. <a href="index.htm">Try again</a>.</p>
</body>
</html>

<?php } ?>