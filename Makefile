CXXFLAGS = -msse4 -O3 -g
all:
	@cd src && ${MAKE} CXXFLAGS="${CXXFLAGS}"

prefix = /usr/local
exec_prefix = ${prefix}
bindir = ${exec_prefix}/bin
install: all
	mkdir -p ${bindir}
	cp bin/tantan ${bindir}

clean:
	@cd src && ${MAKE} clean

tag:
	git tag -m "" `git rev-list HEAD^ | grep -c .`
