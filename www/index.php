
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';
$rforgepkgs = "http://r-forge.r-project.org/R/?group_id=2062";
$tracker = "http://r-forge.r-project.org/tracker/?group_id=2062";
$scm = "http://r-forge.r-project.org/scm/viewvc.php/?root=is2";

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<meta name="keywords" content="stochastic dynamics, inference, hidden Markov model, state space model" />	
	<title><?php echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="80%" cellspacing="0" cellpadding="0">
<tr>
<td>
<a href="http://R-Forge.R-Project.org"><img src="<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> 
</td> 
</tr>
</table>

<table width="80%" align="left" border="0" cellspacing="0" cellpadding="0">
<tr>
<td>
<table border="0" cellspacing="0" cellpadding="0">
<tr>
<td align="right" width="60%">
<h2><span class="emph">is2</span>:<br> 
statistical inference for<br>
<span class="emph">i</span>terated-<span class="emph">s</span>moothing&nbsp;<span class="emph">2</span>nd&nbsp;order</h2>
</td>
<td align="left" width="40%">
<ul>
<li><a href="./index.php?nav=about">About <i>is2</i></a></li>
<li><a href="<?php echo $rforgepkgs;?>">Development Version (Rev. <?php print $svnrevision; ?>)</a></li>
<li><a href="http://lists.r-forge.r-project.org/pipermail/is2-announce/"><i>is2-announce</i> mailing list archives</a></li>
<li><a href="./index.php?nav=news">Package NEWS</a></li>
<li><a href="<?php echo $tracker;?>">Bug reports, feature &amp; support requests</a></li>
<li><a href="<?php echo $scm;?>">Source code repository</a></li>
<li><a href="./index.php?nav=bib">References to the literature</a></li>
<li><a href="./index.php?nav=authors">Authors' homepages</a></li>
</ul>
</td>
</tr>
</table>
</td>
</tr>
<tr>
<td>
<table border="0" cellspacing="0" cellpadding="0">
<tr>
<td>
<?php
$nav = $_REQUEST["nav"];
switch ($nav) {
    case "news":
        $dfile = "content/NEWS.html";
        break;
    case "bib":
        $dfile = "content/refs.htm";
        break;
    case "authors":
        $dfile = "content/authors.htm";
        break;
    default:
        $dfile = "content/about.htm";
        break;
}
include($dfile);
?>
</td>
</tr>
</table>
</td>
</tr>
</table>
</body>
</html>
