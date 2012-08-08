CXXFLAGS = -O3
all:
	@cd src && ${MAKE} CXXFLAGS="${CXXFLAGS}"

prefix = /usr/local
exec_prefix = ${prefix}
bindir = ${exec_prefix}/bin
install: all
	mkdir -p ${bindir}
	cp src/tantan ${bindir}

clean:
	@cd src && ${MAKE} clean

README.html: README.txt
	rst2html README.txt > README.html

log:
	hg log --style changelog > ChangeLog.txt

distdir = tantan-`hg id -n`
dist: README.html log
	@cd src && ${MAKE} version.hh
	rsync -rC --exclude tantan src test Makefile *.txt *.html ${distdir}
	zip -qrm ${distdir} ${distdir}
