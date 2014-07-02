<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:import href="/usr/share/sgml/docbook/xsl-ns-stylesheets/html/chunk.xsl"/>


<xsl:template name="user.footer.content">
  <HR/><TABLE WIDTH="100%"><TR>
  <!--<TD class="copyright">&#x00A9; 2013-2014 Thomas Mejer Hansen.</TD>-->
  <TD class="copyright" ALIGN="right" VALIGN="middlef">This site is hosted by
<A href="http://sourceforge.net/projects/sippi/">sourceforge.net</A></TD>
  </TR></TABLE>
</xsl:template>

<!-- TIC after PREFACE-->
<!--
<xsl:param name="generate.toc"/>
<xsl:param name="process.empty.source.toc" select="1"/>
-->

<xsl:param name="generate.toc" select="'book toc'"/> <!-- ONLY MAIN TOC-->

<xsl:param name="chunk.section.depth" select="2"></xsl:param>
<xsl:param name="chunk.tocs.and.lots" select="0"></xsl:param>
<xsl:param name="chunk.separate.lots" select="0"></xsl:param>

<xsl:param name="use.id.as.filename" select="1"/>

<xsl:param name="xref.with.number.and.title" select="1"/>
<xsl:param name="insert.xref.page.number">yes</xsl:param>

<l:i18n xmlns:l="http://docbook.sourceforge.net/xmlns/l10n/1.0"> 
  <l:l10n language="en"> 
    <l:context name="xref-number-and-title"> 
      <l:template name="appendix" text="Appendix %n: &#8220;%t&#8221;"/> 
      <l:template name="chapter" text="Chapter %n: &#8220;%t&#8221;"/> 
      <l:template name="sect1" text="&#8220;%t&#8221;"/>
    </l:context>    
  </l:l10n>
</l:i18n>

<xsl:param name="html.stylesheet">style_offline.css</xsl:param>

</xsl:stylesheet>

