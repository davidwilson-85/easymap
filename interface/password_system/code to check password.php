<?php //Code for dealing with password====================================================
session_start();
if (isset($_POST['password'])) {
	$_SESSION['password'] = $_POST['password'];
	$password = $_SESSION['password'];
} else {
	@$password = $_SESSION['password'];
}
?>

<?php if ($password == 'elche') { ?>

<<<<<<<<<< HTML CODE TO DISPLAY IF PASSWORD IS CORRECT >>>>>>>>>>

<?php } else { ?>

<html>
<head>
<title></title>
</head>
<body>
Wrong password. <a href="index.htm">Try again</a>
</body>
</html>

<?php } ?>