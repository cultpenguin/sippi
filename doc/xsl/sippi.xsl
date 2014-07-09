<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:import href="/usr/share/xml/docbook/stylesheet/docbook-xsl/html/chunk.xsl"/>


<xsl:template name="user.footer.content">
  <HR/><TABLE WIDTH="100%"><TR>
  <!--<TD class="copyright">&#x00A9; 2013-2014 Thomas Mejer Hansen.</TD>-->
  <TD class="copyright" ALIGN="right" VALIGN="middlef">This site is hosted by
<A href="http://sourceforge.net/projects/sippi/">
sourceforge.net</A></TD>
  </TR></TABLE>
<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-36177321-1', 'sourceforge.net');
  ga('send', 'pageview');

</script>
</xsl:template>

<xsl:param name="chunk.section.depth" select="2"></xsl:param>

<xsl:param name="use.id.as.filename" select="1"/>

<xsl:param name="xref.with.number.and.title" select="1"/>
<xsl:param name="insert.xref.page.number">yes</xsl:param>

<xsl:param name="html.stylesheet">style.css</xsl:param>

</xsl:stylesheet>

