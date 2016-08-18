# Make file for validating/converting sippi.xml
#
#

VERSION=1.5


validate:

	date +"%Y" > yearnow.xml	
	date "+%b %d, %Y" > pubdate.xml	
	echo "${VERSION}" > version.xml
	xmllint --valid --noout sippi.xml
 
dblatex:
	dblatex -T simple -p xsl/dblatex.xsl -b xetex sippi.xml

dblatex_ubuntu:
	dblatex -p xsl/dblatex-ubuntu.xsl -b xetex -o sippi.pdf sippi.xml

pandoc:
	find ./ -name '*.xml' -exec sh -c 'pandoc --from docbook --to markdown -i $1 -o $1.md' _ {} \;

xalan:
	java -jar /usr/share/java/xalan2-2.7.1.jar -in sippi.xml -xsl xsl/sippi.xsl -out sippi.html

saxon:
	saxon-xslt sippi.xml sippi.xsl

fo:
	xsltproc -o sippi.fo /usr/share/xml/docbook/stylesheet/nwalsh/fo/docbook.xsl sippi.xml

pdf:fo
	fop -pdf sippi.pdf -fo sippi.fo

rtf:fo
	fop -rtf sippi.rtf -fo sippi.fo

html:
	make validate
	mkdir -p htmldoc
	mkdir -p htmldoc/figures
	cp figures/*.png htmldoc/figures/.
	cp style*.css htmldoc/.
	xsltproc --xinclude --stringparam html.stylesheet style.css -o htmldoc xsl/sippi.xsl sippi.xml
	mv *.html htmldoc/.

html_cygwin:
	mkdir -p htmldoc
	mkdir -p htmldoc/figures
	cp -r figures htmldoc/.
	cp style*.css htmldoc/.
	xsltproc --xinclude --stringparam -o htmldoc xsl/sippi_cygwin.xsl sippi.xml
	mv *.html htmldoc/.

tidy:
	tidy  -f tidyerrrors.txt  -m  -quiet  htmldoc/*.html


release:
	./make_sippi_release.sh ${VERSION}

clean:
	rm -fr SIPPI
	rm -fr SIPPI_${VERSION}.zip
	rm -f *.htmldoc
	rm -fr htmldoc
	rm -f *.op
	rm -f *.pdf

update_sf:
	./updateSF

all:
	make clean 
	make html 
	make dblatex_ubuntu 
	./updateSF

help:
        # USAGE
        # make linkbook: Use booklayout 
