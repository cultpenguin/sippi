# Make file for validating/converting sippi.xml
#
#

validate:

	date +"%Y" > yearnow.xml	
	date "+%b %d, %Y" > pubdate.xml	
	xmllint --valid --noout sippi.xml

xalan:
	java -jar /usr/share/java/xalan2-2.7.1.jar -in sippi.xml -xsl sippi.xsl -out sippi.html
dblatex:
	dblatex -P doc.collab.show=0 -P latex.output.revhistory=0 -P doc.toc.show=1 sippi.xml

saxon:
	saxon-xslt sippi.xml sippi.xsl

pdf_old:
	make validate
	xsltproc  --output sippi.fo \
	--stringparam fop.extensions 1 \
	file:///usr/share/xml/docbook/stylesheet/nwalsh/fo/docbook.xsl sippi.xml
	fop  -fo  sippi.fo  -pdf  sippi.pdf

pdf_ldp:
	fop -xml sippi.xml -xsl /usr/share/xml/docbook/stylesheet/ldp/fo/tldp-print.xsl  -pdf sippi_ldp.pdf

fo:
	xsltproc -o sippi.fo /usr/share/xml/docbook/stylesheet/nwalsh/fo/docbook.xsl sippi.xml


pdf:fo
	fop -pdf sippi.pdf -fo sippi.fo

rtf:fo
	fop -rtf sippi.rtf -fo sippi.fo

cygwin:	
	xsltproc --xinclude --stringparam html.stylesheet style.css -o htmldoc sippi.xsl sippi.xml

html:
	make validate
	mkdir -p htmldoc
	mkdir -p htmldoc/figures
	cp figures/*.png htmldoc/figures/.
	cp style.css htmldoc/.
	xsltproc --xinclude --stringparam html.stylesheet style.css -o htmldoc sippi.xsl sippi.xml
	mv *.html htmldoc/.

html_cygwin:
	make validate
	mkdir -p htmldoc
	mkdir -p htmldoc/figures
	cp -r figures htmldoc/.
	cp style_offline.css htmldoc/.
	xsltproc --xinclude --stringparam html.stylesheet style.css -o htmldoc sippi_cygwin.xsl sippi.xml
	mv *.html htmldoc/.

	
	

linkbook:
	ln -s sippiBook.xml sippi.xml
linkarticle:
	ln -s sippiArticle.xml sippi.xml
validateall:
	xmllint --valid --noout sippiBook.xml
	xmllint --valid --noout sippiArticle.xml




htmlall:
	xsltproc  --output \
	sippiArticle.html\
	file:///usr/share/xml/docbook/stylesheet/nwalsh/xhtml/docbook.xsl sippiArticle.xml

	xsltproc  --output \
	sippiBook.html\
	file:///usr/share/xml/docbook/stylesheet/nwalsh/xhtml/docbook.xsl sippiBook.xml

allpdf:
	xsltproc  --output sippiBook.fo \
	--stringparam fop.extensions 1 \
	file:///usr/share/xml/docbook/stylesheet/nwalsh/fo/docbook.xsl sippiBook.xml

	fop  -fo  sippiBook.fo  -pdf  sippiBook.pdf

	xsltproc  --output sippiArticle.fo \
	--stringparam fop.extensions 1 \
	file:///usr/share/xml/docbook/stylesheet/nwalsh/fo/docbook.xsl sippiArticle.xml

	fop  -fo  sippiArticle.fo  -pdf  sippiArticle.pdf

clean:
	rm -f *.htmldoc
	rm -fr htmldoc
	rm -f *.op
	rm -f *.pdf

all:html dblatex



help:
	# USAGE
	# make linkbook: Use booklayout 
