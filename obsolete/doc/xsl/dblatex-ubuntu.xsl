<?xml version='1.0' encoding="iso-8859-1"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version='1.0'>

<!-- Let's import our own XSL to override the default behaviour. -->
  <xsl:import href="dblatex.xsl"/>


  <xsl:param name="xetex.font">
    <xsl:text>\setmainfont{Ubuntu}
    </xsl:text>
    <xsl:text>\setsansfont{Ubuntu}
    </xsl:text>
    <xsl:text>\setmonofont{Ubuntu Mono}
    </xsl:text>
  </xsl:param>

 </xsl:stylesheet>
